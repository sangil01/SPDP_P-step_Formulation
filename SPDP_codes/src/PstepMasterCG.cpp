#include "PstepMasterCG.h"
#include "PstepValidInequality.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace spdp {
namespace {

constexpr double kTolerance = 1e-9;

bool double_equal(double lhs, double rhs) {
    return std::fabs(lhs - rhs) <= kTolerance;
}

int to_gurobi_method(GurobiLPMethod method) {
    return static_cast<int>(method);
}

std::string format_double(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << value;
    return out.str();
}

State canonicalize_state(State state) {
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

std::map<NodeId, std::vector<State>> collect_sigma_by_node_from_graph(
    const MultiDiGraph& graph
) {
    std::map<NodeId, std::set<State>> sigma_sets;
    const NodeId end_node_id = graph.end_node_id();

    for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        sigma_sets.emplace(node_id, std::set<State>{});
    }

    for (const EdgeRecord& edge : graph.edges()) {
        if (graph.is_physical_service_node(edge.u)) {
            sigma_sets[edge.u].insert(canonicalize_state(edge.data.start_state));
        }
        if (graph.is_physical_service_node(edge.v)) {
            sigma_sets[edge.v].insert(canonicalize_state(edge.data.end_state));
        }
    }

    std::map<NodeId, std::vector<State>> sigma_by_node;
    for (const auto& entry : sigma_sets) {
        sigma_by_node.emplace(entry.first, std::vector<State>(entry.second.begin(), entry.second.end()));
    }
    return sigma_by_node;
}

void add_phase_one_vi_row(
    CGMasterProblem& problem,
    GRBLinExpr expr,
    double rhs,
    const std::string& name
) {
    GRBVar artificial_var = problem.model->addVar(
        0.0,
        GRB_INFINITY,
        1.0,
        GRB_CONTINUOUS,
        "z_" + name
    );
    problem.artificial_vi_vars.push_back(artificial_var);
    expr += artificial_var;
    problem.model->addConstr(expr >= rhs, name);
}

void add_phase_one_vi_upper_row(
    CGMasterProblem& problem,
    GRBLinExpr expr,
    double rhs,
    const std::string& name
) {
    GRBVar artificial_var = problem.model->addVar(
        0.0,
        GRB_INFINITY,
        1.0,
        GRB_CONTINUOUS,
        "z_" + name
    );
    problem.artificial_vi_vars.push_back(artificial_var);
    expr -= artificial_var;
    problem.model->addConstr(expr <= rhs, name);
}

void add_root_theta_valid_inequalities(
    const SPDPData& data,
    const MultiDiGraph& graph,
    CGMasterProblem& problem,
    const CGMasterOptions& options
) {
    const PstepValidInequalityOptions vi_options{
        options.add_vi_35,
        options.add_vi_36_combined,
        options.vi_36_subset_max_size,
        options.add_vi_request_block_sec,
        options.vi_request_block_sec_max_size,
        options.add_vi_44,
    };
    const std::vector<PstepValidInequalityRow> vi_rows =
        build_pstep_valid_inequality_rows(data, graph, vi_options);
    for (const PstepValidInequalityRow& row : vi_rows) {
        GRBLinExpr expr = 0.0;
        for (const auto& term : row.edge_terms) {
            expr += term.second * problem.theta_vars[term.first];
        }
        if (row.sense == PstepValidInequalitySense::GreaterEqual) {
            add_phase_one_vi_row(problem, expr, row.rhs, row.name);
        } else {
            add_phase_one_vi_upper_row(problem, expr, row.rhs, row.name);
        }
    }
}

std::string build_column_key(
    const std::vector<int>& edge_ids,
    double tau
) {
    std::ostringstream out;
    for (std::size_t idx = 0; idx < edge_ids.size(); ++idx) {
        if (idx > 0U) {
            out << ',';
        }
        out << edge_ids[idx];
    }
    out << '|';
    out << std::fixed << std::setprecision(9) << tau;
    return out.str();
}

}  // namespace

// x_r이 없는 상태에서 artificial variable를 없애기 위한 phase 1 master problem을 생성한다.
CGMasterProblem build_cg_master_problem(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const CGMasterOptions& options
) {
    CGMasterProblem problem;
    problem.phase = CGPhase::PhaseI;

    problem.env = std::make_unique<GRBEnv>(true);
    if (!options.gurobi_log_path.empty()) {
        problem.env->set(GRB_StringParam_LogFile, options.gurobi_log_path);
    }
    problem.env->start();

    problem.model = std::make_unique<GRBModel>(*problem.env);
    problem.model->set(GRB_StringAttr_ModelName, "node_cg_master");
    problem.model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    problem.model->set(GRB_IntParam_Method, to_gurobi_method(options.initial_lp_method));
    if (options.gurobi_threads >= 0) {
        problem.model->set(GRB_IntParam_Threads, options.gurobi_threads);
    }
    if (options.solver_time_limit >= 0.0) {
        problem.model->set(GRB_DoubleParam_TimeLimit, options.solver_time_limit);
    }

    const NodeId end_node_id = graph.end_node_id();
    problem.sigma_by_node = collect_sigma_by_node_from_graph(graph);
    problem.fixed_theta_value_by_edge.assign(graph.number_of_edges(), -1);
    for (std::size_t edge_id = 0;
         edge_id < options.fixed_theta_value_by_edge.size() &&
         edge_id < graph.number_of_edges();
         ++edge_id) {
        problem.fixed_theta_value_by_edge[edge_id] = options.fixed_theta_value_by_edge[edge_id];
    }

    // theta_e variables: node LP relaxation에서는 [0,1] continuous.
    problem.theta_vars.reserve(graph.number_of_edges());
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        double lower_bound = 0.0;
        double upper_bound = 1.0;
        if (problem.fixed_theta_value_by_edge[edge_id] == 0) {
            lower_bound = 0.0;
            upper_bound = 0.0;
        } else if (problem.fixed_theta_value_by_edge[edge_id] == 1) {
            lower_bound = 1.0;
            upper_bound = 1.0;
        }
        problem.theta_vars.push_back(
            problem.model->addVar(
                lower_bound,
                upper_bound,
                0.0,
                GRB_CONTINUOUS,
                "theta_" + std::to_string(edge_id)
            )
        );
    }

    // Phase I visit artificial z_i variables.
    for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        problem.artificial_visit_vars.push_back(
            problem.model->addVar(
                0.0,
                GRB_INFINITY,
                1.0,
                GRB_CONTINUOUS,
                "z_visit_" + std::to_string(node_id)
            )
        );
    }
    problem.model->update();

    int artificial_index = 0;
    for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }
        GRBLinExpr visit_expr = 0.0;
        visit_expr += problem.artificial_visit_vars[static_cast<std::size_t>(artificial_index)];
        problem.visit_rows.emplace(
            node_id,
            problem.model->addConstr(
                visit_expr == 2.0,
                "visit_" + std::to_string(node_id)
            )
        );
        ++artificial_index;
    }

    for (const auto& sigma_entry : problem.sigma_by_node) {
        const NodeId node_id = sigma_entry.first;
        const std::vector<State>& sigma_values = sigma_entry.second;

        for (std::size_t sigma_idx = 0; sigma_idx < sigma_values.size(); ++sigma_idx) {
            const NodeStateKey key{node_id, canonicalize_state(sigma_values[sigma_idx])};

            const std::string suffix =
                std::to_string(node_id) + "_" + std::to_string(sigma_idx);

            problem.artificial_state_pos_vars.emplace(
                key,
                problem.model->addVar(
                    0.0,
                    GRB_INFINITY,
                    1.0,
                    GRB_CONTINUOUS,
                    "z_state_pos_" + suffix
                )
            );
            problem.artificial_state_neg_vars.emplace(
                key,
                problem.model->addVar(
                    0.0,
                    GRB_INFINITY,
                    1.0,
                    GRB_CONTINUOUS,
                    "z_state_neg_" + suffix
                )
            );
            problem.artificial_time_vars.emplace(
                key,
                problem.model->addVar(
                    0.0,
                    GRB_INFINITY,
                    1.0,
                    GRB_CONTINUOUS,
                    "h_time_" + suffix
                )
            );

            GRBLinExpr state_expr = 0.0;
            state_expr += problem.artificial_state_pos_vars.at(key);
            state_expr -= problem.artificial_state_neg_vars.at(key);
            problem.state_rows.emplace(
                key,
                problem.model->addConstr(
                    state_expr == 0.0,
                    "state_" + suffix
                )
            );

            GRBLinExpr time_expr = 0.0;
            time_expr += problem.artificial_time_vars.at(key);
            problem.time_rows.emplace(
                key,
                problem.model->addConstr(
                    time_expr >= 0.0,
                    "time_" + suffix
                )
            );
        }
    }

    problem.edge_rows.reserve(graph.number_of_edges());
    for (std::size_t edge_id = 0; edge_id < graph.number_of_edges(); ++edge_id) {
        GRBLinExpr edge_expr = -problem.theta_vars[edge_id];
        if (problem.fixed_theta_value_by_edge[edge_id] == 1) {
            GRBVar artificial_var = problem.model->addVar(
                0.0,
                GRB_INFINITY,
                1.0,
                GRB_CONTINUOUS,
                "z_theta_fix_" + std::to_string(edge_id)
            );
            problem.artificial_fixed_theta_vars.emplace(edge_id, artificial_var);
            edge_expr += artificial_var;
        }
        problem.edge_rows.push_back(
            problem.model->addConstr(
                edge_expr == 0.0,
                "edge_" + std::to_string(edge_id)
            )
        );
    }

    add_root_theta_valid_inequalities(data, graph, problem, options);

    return problem;
}

std::size_t add_columns_to_cg_master(
    CGMasterProblem& problem,
    const std::vector<CGColumn>& columns
) {
    std::size_t added_count = 0;

    for (const CGColumn& column : columns) {
        const std::string column_key = build_column_key(column.edge_ids, column.tau);
        if (problem.column_id_by_key.find(column_key) != problem.column_id_by_key.end()) {
            continue;
        }

        GRBColumn grb_column;

        for (const auto& entry : column.visit_coefficients) {
            const auto found = problem.visit_rows.find(entry.first);
            if (found != problem.visit_rows.end()) {
                grb_column.addTerm(static_cast<double>(entry.second), found->second);
            }
        }

        for (const auto& entry : column.state_coefficients) {
            const auto found = problem.state_rows.find(entry.first);
            if (found != problem.state_rows.end()) {
                grb_column.addTerm(static_cast<double>(entry.second), found->second);
            }
        }

        for (const auto& entry : column.time_coefficients) {
            const auto found = problem.time_rows.find(entry.first);
            if (found != problem.time_rows.end()) {
                grb_column.addTerm(entry.second, found->second);
            }
        }

        for (int edge_id : column.edge_incidence) {
            grb_column.addTerm(1.0, problem.edge_rows[static_cast<std::size_t>(edge_id)]);
        }

        const double objective_coefficient =
            (problem.phase == CGPhase::PhaseII) ? column.total_cost : 0.0;
        const int column_id = static_cast<int>(problem.columns.size());
        GRBVar x_var = problem.model->addVar(
            0.0,
            GRB_INFINITY,
            objective_coefficient,
            GRB_CONTINUOUS,
            grb_column,
            "x_" + std::to_string(column_id)
        );

        CGColumn stored_column = column;
        stored_column.id = column_id;
        problem.x_vars.push_back(x_var);
        problem.columns.push_back(std::move(stored_column));
        problem.column_id_by_key.emplace(column_key, column_id);
        ++added_count;
    }

    if (added_count > 0U) {
        problem.model->update();
    }

    return added_count;
}

void solve_cg_master_lp(CGMasterProblem& problem) {
    solve_cg_master_lp(problem, GurobiLPMethod::Automatic);
}

void solve_cg_master_lp(
    CGMasterProblem& problem,
    GurobiLPMethod method
) {
    if (problem.model == nullptr) {
        throw std::runtime_error("CG master model must exist before solving.");
    }
    if (method != GurobiLPMethod::Automatic) {
        problem.model->set(GRB_IntParam_Method, to_gurobi_method(method));
    }
    problem.model->optimize();
}

CGDualSolution extract_cg_master_duals(const CGMasterProblem& problem) {
    CGDualSolution duals;

    for (const auto& entry : problem.visit_rows) {
        duals.visit_duals.emplace(entry.first, entry.second.get(GRB_DoubleAttr_Pi));
    }

    for (const auto& entry : problem.state_rows) {
        duals.state_duals.emplace(entry.first, entry.second.get(GRB_DoubleAttr_Pi));
    }

    for (const auto& entry : problem.time_rows) {
        duals.time_duals.emplace(entry.first, entry.second.get(GRB_DoubleAttr_Pi));
    }

    duals.edge_duals.reserve(problem.edge_rows.size());
    for (const GRBConstr& constr : problem.edge_rows) {
        duals.edge_duals.push_back(constr.get(GRB_DoubleAttr_Pi));
    }

    return duals;
}

CGLPSnapshot capture_cg_master_snapshot(const CGMasterProblem& problem) {
    CGLPSnapshot snapshot;
    snapshot.gurobi_status = problem.model->get(GRB_IntAttr_Status); 
    snapshot.runtime_seconds = problem.model->get(GRB_DoubleAttr_Runtime);

    const int solution_count = problem.model->get(GRB_IntAttr_SolCount);
    snapshot.has_primal_solution = solution_count > 0;
    if (solution_count > 0) {
        try {
            snapshot.objective_value = problem.model->get(GRB_DoubleAttr_ObjVal);
        } catch (const GRBException&) {
        }
        try {
            snapshot.objective_bound = problem.model->get(GRB_DoubleAttr_ObjBound);
        } catch (const GRBException&) {
        }

        for (const GRBVar& artificial_var : problem.artificial_visit_vars) {
            snapshot.artificial_sum += artificial_var.get(GRB_DoubleAttr_X);
        }
        for (const auto& entry : problem.artificial_state_pos_vars) {
            snapshot.artificial_sum += entry.second.get(GRB_DoubleAttr_X);
        }
        for (const auto& entry : problem.artificial_state_neg_vars) {
            snapshot.artificial_sum += entry.second.get(GRB_DoubleAttr_X);
        }
        for (const auto& entry : problem.artificial_time_vars) {
            snapshot.artificial_sum += entry.second.get(GRB_DoubleAttr_X);
        }
        for (const auto& entry : problem.artificial_fixed_theta_vars) {
            snapshot.artificial_sum += entry.second.get(GRB_DoubleAttr_X);
        }
        for (const GRBVar& artificial_var : problem.artificial_vi_vars) {
            snapshot.artificial_sum += artificial_var.get(GRB_DoubleAttr_X);
        }
        for (const GRBVar& x_var : problem.x_vars) {
            if (x_var.get(GRB_DoubleAttr_X) > 1e-6) {
                ++snapshot.positive_x_count;
            }
        }
        for (const GRBVar& theta_var : problem.theta_vars) {
            if (theta_var.get(GRB_DoubleAttr_X) > 1e-6) {
                ++snapshot.active_theta_count;
            }
        }
    } else {
        try {
            snapshot.objective_bound = problem.model->get(GRB_DoubleAttr_ObjBound);
        } catch (const GRBException&) {
        }
    }

    return snapshot;
}

void switch_cg_master_to_phase_two(CGMasterProblem& problem) {
    problem.phase = CGPhase::PhaseII;

    for (GRBVar& artificial_var : problem.artificial_visit_vars) {
        artificial_var.set(GRB_DoubleAttr_Obj, 0.0);
        artificial_var.set(GRB_DoubleAttr_UB, 0.0);
    }
    for (auto& entry : problem.artificial_state_pos_vars) {
        entry.second.set(GRB_DoubleAttr_Obj, 0.0);
        entry.second.set(GRB_DoubleAttr_UB, 0.0);
    }
    for (auto& entry : problem.artificial_state_neg_vars) {
        entry.second.set(GRB_DoubleAttr_Obj, 0.0);
        entry.second.set(GRB_DoubleAttr_UB, 0.0);
    }
    for (auto& entry : problem.artificial_time_vars) {
        entry.second.set(GRB_DoubleAttr_Obj, 0.0);
        entry.second.set(GRB_DoubleAttr_UB, 0.0);
    }
    for (auto& entry : problem.artificial_fixed_theta_vars) {
        entry.second.set(GRB_DoubleAttr_Obj, 0.0);
        entry.second.set(GRB_DoubleAttr_UB, 0.0);
    }
    for (GRBVar& artificial_var : problem.artificial_vi_vars) {
        artificial_var.set(GRB_DoubleAttr_Obj, 0.0);
        artificial_var.set(GRB_DoubleAttr_UB, 0.0);
    }
    for (std::size_t column_id = 0; column_id < problem.x_vars.size(); ++column_id) {
        problem.x_vars[column_id].set(
            GRB_DoubleAttr_Obj,
            problem.columns[column_id].total_cost
        );
    }

    problem.model->update();
}

void write_cg_master_snapshot(
    std::ostream& out,
    const CGLPSnapshot& snapshot
) {
    out << "[cg-master] status=" << snapshot.gurobi_status
        << " has_primal=" << (snapshot.has_primal_solution ? 1 : 0)
        << " obj=" << format_double(snapshot.objective_value)
        << " bound=" << format_double(snapshot.objective_bound)
        << " runtime=" << format_double(snapshot.runtime_seconds)
        << " artificial_sum=" << format_double(snapshot.artificial_sum)
        << " positive_x=" << snapshot.positive_x_count
        << " active_theta=" << snapshot.active_theta_count
        << '\n';
}

}  // namespace spdp
