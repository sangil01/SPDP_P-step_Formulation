#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "GenMultiGraph.h"
#include "PstepBnP.h"
#include "PstepFormulation.h"
#include "ReadData.h"

namespace {

struct CliOptions {
    std::string instance = "A0.dat";
    int p = 2;
    double solver_time_limit = 3600.0;
    int gurobi_threads = -1;
    std::string solver_mode = "enumeration";
    std::string node_cg_phase_one_mode = "exact-cg";
    std::string vi_formulation = "theta";
    int dump_psteps = 10;
    int validate_psteps = 0;
    int solve_model = 1;
    int prune_infeasible_edges = 1;
    int prune_dominated_edges = 1;
    int prune_symmetry_40 = 1;
    int prune_symmetry_41 = 1;
    int prune_pickup_symmetry_43 = 1;
    int prune_delivery_symmetry_43 = 1;
    int add_vi_35 = 0;
    int add_vi_36 = 0;
    int add_vi_44 = 0;
    int cg_max_iterations_per_phase = 1000;
    int cg_max_columns_per_start = 16;
    int cg_max_total_columns_per_round = 1024;
    double cg_reduced_cost_tolerance = -1e-6;
};

void print_usage(const char* executable) {
    std::cerr << "Usage: " << executable
              << " [instance] [--p N] [--solver-time-limit T] [--gurobi-threads N]"
              << " [--solver-mode enumeration|root-cg]"
              << " [--node-cg-phase1-mode exact-cg|heuristic-cg]"
              << " [--vi-formulation theta|x] [--dump-psteps N] [--validate-psteps 0|1] [--solve 0|1]"
              << " [--prune-infeasible-edges 0|1] [--prune-dominated-edges 0|1]"
              << " [--prune-symmetry-40 0|1] [--prune-symmetry-41 0|1]"
              << " [--prune-pickup-symmetry-43 0|1]"
              << " [--prune-delivery-symmetry-43 0|1]"
              << " [--add-vi-35 0|1] [--add-vi-36 0|1]"
              << " [--add-vi-44 0|1]"
              << " [--cg-max-iterations-per-phase N]"
              << " [--cg-max-columns-per-start N]"
              << " [--cg-max-total-columns-per-round N]"
              << " [--cg-reduced-cost-tolerance T]\n";
}

int parse_int(const std::string& value, const std::string& field_name) {
    try {
        std::size_t consumed = 0;
        const int parsed = std::stoi(value, &consumed);
        if (consumed != value.size()) {
            throw std::invalid_argument("Trailing characters");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid integer for " + field_name + ": " + value);
    }
}

double parse_double(const std::string& value, const std::string& field_name) {
    try {
        std::size_t consumed = 0;
        const double parsed = std::stod(value, &consumed);
        if (consumed != value.size()) {
            throw std::invalid_argument("Trailing characters");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid numeric value for " + field_name + ": " + value);
    }
}

std::string parse_vi_formulation(const std::string& value) {
    if (value == "theta" || value == "x") {
        return value;
    }
    throw std::runtime_error("Invalid value for --vi-formulation: " + value + " (expected theta or x)");
}

std::string parse_solver_mode(const std::string& value) {
    if (value == "enumeration" || value == "root-cg") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for --solver-mode: " + value +
        " (expected enumeration or root-cg)"
    );
}

std::string parse_node_cg_phase_one_mode(const std::string& value) {
    if (value == "exact-cg" || value == "heuristic-cg") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for node CG phase-one mode: " + value +
        " (expected exact-cg or heuristic-cg)"
    );
}

spdp::VIFormulation to_vi_formulation(const std::string& value) {
    if (value == "theta") {
        return spdp::VIFormulation::Theta;
    }
    if (value == "x") {
        return spdp::VIFormulation::X;
    }
    throw std::runtime_error("Unsupported VI formulation: " + value);
}

CliOptions parse_cli(int argc, char** argv) {
    CliOptions options;
    bool instance_set = false;

    for (int idx = 1; idx < argc; ++idx) {
        const std::string arg = argv[idx];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            std::exit(0);
        }

        if (arg == "--p") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--p requires a value.");
            }
            options.p = parse_int(argv[++idx], "--p");
            continue;
        }

        if (arg == "--solver-mode") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--solver-mode requires a value.");
            }
            options.solver_mode = parse_solver_mode(argv[++idx]);
            continue;
        }

        if (arg == "--node-cg-phase1-mode" || arg == "--root-phase1-mode") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.node_cg_phase_one_mode = parse_node_cg_phase_one_mode(argv[++idx]);
            continue;
        }

        if (arg == "--dump-psteps") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--dump-psteps requires a value.");
            }
            options.dump_psteps = parse_int(argv[++idx], "--dump-psteps");
            continue;
        }

        if (arg == "--vi-formulation") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--vi-formulation requires a value.");
            }
            options.vi_formulation = parse_vi_formulation(argv[++idx]);
            continue;
        }

        if (arg == "--validate-psteps") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--validate-psteps requires a value.");
            }
            options.validate_psteps = parse_int(argv[++idx], "--validate-psteps");
            if (options.validate_psteps != 0 && options.validate_psteps != 1) {
                throw std::runtime_error("--validate-psteps must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--solver-time-limit") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--solver-time-limit requires a value.");
            }
            options.solver_time_limit = parse_double(argv[++idx], "--solver-time-limit");
            continue;
        }

        if (arg == "--cg-max-iterations-per-phase") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--cg-max-iterations-per-phase requires a value.");
            }
            options.cg_max_iterations_per_phase =
                parse_int(argv[++idx], "--cg-max-iterations-per-phase");
            if (options.cg_max_iterations_per_phase <= 0) {
                throw std::runtime_error("--cg-max-iterations-per-phase must be positive.");
            }
            continue;
        }

        if (arg == "--cg-max-columns-per-start") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--cg-max-columns-per-start requires a value.");
            }
            options.cg_max_columns_per_start =
                parse_int(argv[++idx], "--cg-max-columns-per-start");
            if (options.cg_max_columns_per_start <= 0) {
                throw std::runtime_error("--cg-max-columns-per-start must be positive.");
            }
            continue;
        }

        if (arg == "--cg-max-total-columns-per-round") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--cg-max-total-columns-per-round requires a value.");
            }
            options.cg_max_total_columns_per_round =
                parse_int(argv[++idx], "--cg-max-total-columns-per-round");
            if (options.cg_max_total_columns_per_round <= 0) {
                throw std::runtime_error("--cg-max-total-columns-per-round must be positive.");
            }
            continue;
        }

        if (arg == "--cg-reduced-cost-tolerance") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--cg-reduced-cost-tolerance requires a value.");
            }
            options.cg_reduced_cost_tolerance =
                parse_double(argv[++idx], "--cg-reduced-cost-tolerance");
            if (options.cg_reduced_cost_tolerance > 0.0) {
                throw std::runtime_error("--cg-reduced-cost-tolerance must be nonpositive.");
            }
            continue;
        }

        if (arg == "--gurobi-threads") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--gurobi-threads requires a value.");
            }
            options.gurobi_threads = parse_int(argv[++idx], "--gurobi-threads");
            if (options.gurobi_threads < 0) {
                throw std::runtime_error("--gurobi-threads must be nonnegative.");
            }
            continue;
        }

        if (arg == "--solve") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--solve requires a value.");
            }
            options.solve_model = parse_int(argv[++idx], "--solve");
            if (options.solve_model != 0 && options.solve_model != 1) {
                throw std::runtime_error("--solve must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--prune-dominated-edges") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--prune-dominated-edges requires a value.");
            }
            options.prune_dominated_edges = parse_int(argv[++idx], "--prune-dominated-edges");
            if (options.prune_dominated_edges != 0 && options.prune_dominated_edges != 1) {
                throw std::runtime_error("--prune-dominated-edges must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--prune-infeasible-edges") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--prune-infeasible-edges requires a value.");
            }
            options.prune_infeasible_edges = parse_int(argv[++idx], "--prune-infeasible-edges");
            if (options.prune_infeasible_edges != 0 && options.prune_infeasible_edges != 1) {
                throw std::runtime_error("--prune-infeasible-edges must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--prune-symmetry-40") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--prune-symmetry-40 requires a value.");
            }
            options.prune_symmetry_40 = parse_int(argv[++idx], "--prune-symmetry-40");
            if (options.prune_symmetry_40 != 0 && options.prune_symmetry_40 != 1) {
                throw std::runtime_error("--prune-symmetry-40 must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--prune-symmetry-41") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--prune-symmetry-41 requires a value.");
            }
            options.prune_symmetry_41 = parse_int(argv[++idx], "--prune-symmetry-41");
            if (options.prune_symmetry_41 != 0 && options.prune_symmetry_41 != 1) {
                throw std::runtime_error("--prune-symmetry-41 must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--prune-pickup-symmetry-43") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--prune-pickup-symmetry-43 requires a value.");
            }
            options.prune_pickup_symmetry_43 =
                parse_int(argv[++idx], "--prune-pickup-symmetry-43");
            if (options.prune_pickup_symmetry_43 != 0 &&
                options.prune_pickup_symmetry_43 != 1) {
                throw std::runtime_error("--prune-pickup-symmetry-43 must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--prune-delivery-symmetry-43") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--prune-delivery-symmetry-43 requires a value.");
            }
            options.prune_delivery_symmetry_43 =
                parse_int(argv[++idx], "--prune-delivery-symmetry-43");
            if (options.prune_delivery_symmetry_43 != 0 &&
                options.prune_delivery_symmetry_43 != 1) {
                throw std::runtime_error("--prune-delivery-symmetry-43 must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--add-vi-35") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--add-vi-35 requires a value.");
            }
            options.add_vi_35 = parse_int(argv[++idx], "--add-vi-35");
            if (options.add_vi_35 != 0 && options.add_vi_35 != 1) {
                throw std::runtime_error("--add-vi-35 must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--add-vi-36") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--add-vi-36 requires a value.");
            }
            options.add_vi_36 = parse_int(argv[++idx], "--add-vi-36");
            if (options.add_vi_36 != 0 && options.add_vi_36 != 1) {
                throw std::runtime_error("--add-vi-36 must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--add-vi-44") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--add-vi-44 requires a value.");
            }
            options.add_vi_44 = parse_int(argv[++idx], "--add-vi-44");
            if (options.add_vi_44 != 0 && options.add_vi_44 != 1) {
                throw std::runtime_error("--add-vi-44 must be 0 or 1.");
            }
            continue;
        }

        if (!arg.empty() && arg[0] == '-') {
            throw std::runtime_error("Unknown option: " + arg);
        }

        if (!instance_set) {
            options.instance = arg;
            instance_set = true;
        } else {
            throw std::runtime_error("Too many positional arguments.");
        }
    }

    return options;
}

std::string format_double(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

std::string format_gurobi_threads(int threads) {
    if (threads < 0) {
        return "default(auto)";
    }
    return std::to_string(threads);
}

std::string format_sequence_pi(const std::vector<int>& sequence_pi) {
    if (sequence_pi.empty()) {
        return "()";
    }

    std::ostringstream out;
    out << "(";
    for (std::size_t idx = 0; idx < sequence_pi.size(); ++idx) {
        if (idx > 0U) {
            out << ", ";
        }
        out << sequence_pi[idx];
    }
    if (sequence_pi.size() == 1U) {
        out << ",";
    }
    out << ")";
    return out.str();
}

std::filesystem::path project_root_path() {
    const std::filesystem::path source_path(__FILE__);
    return source_path.parent_path().parent_path().parent_path();
}

std::filesystem::path build_log_output_path(const std::string& instance, int p) {
    const std::filesystem::path instance_path(instance);
    const std::string output_name =
        instance_path.stem().string() + "_p" + std::to_string(p) + "_log.txt";
    return project_root_path() / "SPDP_output" / output_name;
}

std::filesystem::path build_solution_output_path(const std::string& instance, int p) {
    const std::filesystem::path instance_path(instance);
    const std::string output_name =
        instance_path.stem().string() + "_p" + std::to_string(p) + "_sol.txt";
    return project_root_path() / "SPDP_output" / output_name;
}

std::filesystem::path build_gurobi_log_path(const std::string& instance, int p) {
    const std::filesystem::path instance_path(instance);
    const std::string output_name =
        instance_path.stem().string() + "_p" + std::to_string(p) + "_gurobi.log";
    return project_root_path() / "SPDP_output" / output_name;
}

void print_instance_summary(
    std::ostream& out,
    const CliOptions& args,
    const spdp::SPDPData& data,
    const spdp::MultiDiGraph& graph
) {
    out << "[main] Loaded instance: " << args.instance << '\n';
    out << "[main] p: " << args.p << '\n';
    out << "[main] Solver mode: " << args.solver_mode << '\n';
    out << "[main] Node CG phase-1 mode: " << args.node_cg_phase_one_mode << '\n';
    out << "[main] Solver time limit: " << format_double(args.solver_time_limit) << '\n';
    out << "[main] Gurobi threads: " << format_gurobi_threads(args.gurobi_threads) << '\n';
    out << "[main] VI formulation: " << args.vi_formulation << '\n';
    out << "[main] Locations: " << data.locations << '\n';
    out << "[main] Fixed vehicle cost: " << format_double(data.fixed_vehicle_cost) << '\n';
    out << "[main] Time pick-up: " << format_double(data.time_pickup) << '\n';
    out << "[main] Time empty: " << format_double(data.time_empty) << '\n';
    out << "[main] Time delivery: " << format_double(data.time_delivery) << '\n';
    out << "[main] Time limit: " << format_double(data.time_limit) << '\n';
    out << "[main] Requests: " << data.requests.size() << '\n';
    out << "[main] Time matrix size: " << data.time.size() << "x"
        << (data.time.empty() ? 0 : data.time.front().size()) << '\n';
    out << "[main] Distance matrix size: " << data.distance.size() << "x"
        << (data.distance.empty() ? 0 : data.distance.front().size()) << '\n';
    out << "[main] Node count: " << graph.number_of_nodes() << '\n';
    out << "[main] Edge count: " << graph.number_of_edges() << '\n';
    out << "[main] prune_infeasible_edges: " << args.prune_infeasible_edges << '\n';
    out << "[main] prune_dominated_edges: " << args.prune_dominated_edges << '\n';
    out << "[main] prune_symmetry_40: " << args.prune_symmetry_40 << '\n';
    out << "[main] prune_symmetry_41: " << args.prune_symmetry_41 << '\n';
    out << "[main] prune_pickup_symmetry_43: " << args.prune_pickup_symmetry_43 << '\n';
    out << "[main] prune_delivery_symmetry_43: " << args.prune_delivery_symmetry_43 << '\n';
    out << "[main] add_vi_35: " << args.add_vi_35 << '\n';
    out << "[main] add_vi_36: " << args.add_vi_36 << '\n';
    out << "[main] add_vi_44: " << args.add_vi_44 << '\n';
    out << "[main] cg_max_iterations_per_phase: " << args.cg_max_iterations_per_phase << '\n';
    out << "[main] cg_max_columns_per_start: " << args.cg_max_columns_per_start << '\n';
    out << "[main] cg_max_total_columns_per_round: " << args.cg_max_total_columns_per_round << '\n';
    out << "[main] cg_reduced_cost_tolerance: "
        << std::scientific << std::setprecision(6) << args.cg_reduced_cost_tolerance << '\n';
    out << std::defaultfloat;
}

void print_selected_edge_info(
    std::ostream& out,
    const spdp::MultiDiGraph& graph,
    bool show_full_edge_list = false
) {
    std::set<std::pair<int, int>> unique_pairs;
    if (show_full_edge_list) {
        for (const spdp::EdgeRecord& edge : graph.edges()) {
            unique_pairs.insert({edge.u, edge.v});
        }
    }

    out << "[main] Full edge list for selected node pairs:\n";
    for (const auto& pair : unique_pairs) {
        const std::vector<spdp::EdgeRecord> pair_edges = graph.edges_between(pair.first, pair.second);
        out << "  Pair (" << pair.first << ", " << pair.second << ") -> " << pair_edges.size()
            << " edges\n";

        for (const spdp::EdgeRecord& edge : pair_edges) {
            out << "    key=" << edge.key
                << " | pi=" << format_sequence_pi(edge.data.sequence_pi)
                << " | time=" << format_double(edge.data.time)
                << " | cost=" << format_double(edge.data.cost)
                << " | start=" << spdp::state_to_str(edge.data.start_state)
                << " | end=" << spdp::state_to_str(edge.data.end_state)
                << '\n';
        }
    }
}

void print_pstep_summary(
    std::ostream& out,
    const spdp::CompactPStepArtifacts& artifacts
) {
    out << "[main] Raw feasible p-step count: " << artifacts.raw_paths.size() << '\n';
    out << "[main] Compact p-step count: " << artifacts.compact_psteps.size() << '\n';
    out << "[main] Sigma_i state counts:\n";
    for (const auto& entry : artifacts.coefficients.sigma_by_node) {
        out << "  node " << entry.first << " -> " << entry.second.size() << '\n';
    }
}

void print_solution_summary(
    std::ostream& out,
    const spdp::CompactMasterProblem& problem
) {
    const int status = problem.model->get(GRB_IntAttr_Status);
    out << "[main] Gurobi status: " << status << '\n';

    const double lower_bound = problem.model->get(GRB_DoubleAttr_ObjBound);

    const int solution_count = problem.model->get(GRB_IntAttr_SolCount);
    if (solution_count <= 0) {
        out << "[main] No feasible solution available.\n";
        out << "[main] Lower bound: " << format_double(lower_bound) << '\n';
        return;
    }

    const double objective_value = problem.model->get(GRB_DoubleAttr_ObjVal);
    const double gap_percent =
        std::abs(objective_value) > 1e-9
            ? 100.0 * std::abs(objective_value - lower_bound) / std::abs(objective_value)
            : 0.0;

    out << "[main] Objective value: " << format_double(objective_value) << '\n';
    out << "[main] Lower bound: " << format_double(lower_bound) << '\n';
    out << "[main] Gap (%): " << format_double(gap_percent) << '\n';
    out << "[main] Solver runtime (sec): "
        << format_double(problem.model->get(GRB_DoubleAttr_Runtime)) << '\n';

    int positive_x_count = 0;
    for (const GRBVar& var : problem.x_vars) {
        if (var.get(GRB_DoubleAttr_X) > 1e-6) {
            ++positive_x_count;
        }
    }

    int active_theta_count = 0;
    for (const GRBVar& var : problem.theta_vars) {
        if (var.get(GRB_DoubleAttr_X) > 0.5) {
            ++active_theta_count;
        }
    }

    out << "[main] Positive x_r count: " << positive_x_count << '\n';
    out << "[main] Active theta_e count: " << active_theta_count << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const CliOptions args = parse_cli(argc, argv);
        const spdp::SPDPData data = spdp::read_spdp_data(args.instance);

        const std::filesystem::path log_output_path = build_log_output_path(args.instance, args.p);
        const std::filesystem::path solution_output_path =
            build_solution_output_path(args.instance, args.p);
        const std::filesystem::path gurobi_log_path =
            build_gurobi_log_path(args.instance, args.p);
        std::filesystem::create_directories(log_output_path.parent_path());

        std::ofstream output_file(log_output_path);
        if (!output_file) {
            throw std::runtime_error("Failed to open output file: " + log_output_path.string());
        }

        std::ofstream solution_file(solution_output_path);
        if (!solution_file) {
            throw std::runtime_error(
                "Failed to open solution output file: " + solution_output_path.string()
            );
        }

        const spdp::GraphBuildOptions graph_build_options{
            args.prune_infeasible_edges == 1,
            args.prune_dominated_edges == 1,
            args.prune_symmetry_40 == 1,
            args.prune_symmetry_41 == 1,
        };
        const spdp::MultiDiGraph graph = spdp::build_multigraph(
            data,
            graph_build_options,
            &output_file
        );

        print_instance_summary(output_file, args, data, graph);
        print_selected_edge_info(output_file, graph);

        std::ofstream gurobi_log_file(gurobi_log_path, std::ios::trunc);
        if (!gurobi_log_file) {
            throw std::runtime_error("Failed to initialize Gurobi log file: " + gurobi_log_path.string());
        }
        gurobi_log_file.close();

        if (args.solver_mode == "root-cg") {
            if (args.add_vi_35 == 1 || args.add_vi_36 == 1 || args.add_vi_44 == 1) {
                throw std::runtime_error(
                    "root-cg mode does not yet support VI-35/36/44. "
                    "The seedless Phase-I RMP becomes infeasible without dedicated artificials "
                    "for those inequalities. Please set them to 0."
                );
            }
            if (args.vi_formulation != "theta") {
                throw std::runtime_error(
                    "root-cg mode currently supports only --vi-formulation theta."
                );
            }

            if (args.solve_model == 1) {
                const spdp::NodeCGOptions cg_options{
                    args.p,
                    data.time_limit,
                    args.solver_time_limit,
                    args.gurobi_threads,
                    static_cast<std::size_t>(args.cg_max_iterations_per_phase),
                    static_cast<std::size_t>(args.cg_max_columns_per_start),
                    static_cast<std::size_t>(args.cg_max_total_columns_per_round),
                    args.cg_reduced_cost_tolerance,
                    args.prune_pickup_symmetry_43 == 1,
                    args.prune_delivery_symmetry_43 == 1,
                    gurobi_log_path.string(),
                    args.node_cg_phase_one_mode == "exact-cg"
                        ? spdp::NodeCGPhaseOneMode::ExactCG
                        : spdp::NodeCGPhaseOneMode::HeuristicCG,
                };
                const spdp::NodeCGResult cg_result =
                    spdp::solve_node_column_generation(data, graph, cg_options, &output_file);
                spdp::write_node_cg_summary(output_file, cg_result);
                spdp::write_node_cg_solution(solution_file, cg_result);
            } else {
                output_file << "[main] Solve skipped by CLI option.\n";
                solution_file << "!! Print 2\n";
                solution_file << "Solution:\n";
                solution_file << "  Solve skipped by CLI option\n";
                solution_file << "Solution done \n";
            }
        } else {
            const spdp::CompactPStepOptions options{
                args.p,
                data.time_limit,
                static_cast<std::size_t>(std::max(args.dump_psteps, 0)),
                args.validate_psteps == 1,
                args.prune_pickup_symmetry_43 == 1,
                args.prune_delivery_symmetry_43 == 1,
            };
            const auto compact_pstep_build_start = std::chrono::steady_clock::now();
            const spdp::CompactPStepArtifacts artifacts =
                spdp::build_compact_pstep_artifacts(graph, options, &output_file);
            const auto compact_pstep_build_end = std::chrono::steady_clock::now();
            const double compact_pstep_build_seconds =
                std::chrono::duration<double>(compact_pstep_build_end - compact_pstep_build_start)
                    .count();

            print_pstep_summary(output_file, artifacts);
            output_file << "[main] Compact p-step set build time (sec): "
                        << format_double(compact_pstep_build_seconds) << '\n';
            spdp::dump_compact_psteps(output_file, artifacts.compact_psteps, options.dump_limit);

            const spdp::CompactMasterBuildOptions master_build_options{
                gurobi_log_path.string(),
                args.add_vi_35 == 1,
                args.add_vi_36 == 1,
                args.add_vi_44 == 1,
                to_vi_formulation(args.vi_formulation),
            };
            spdp::CompactMasterProblem problem = spdp::build_compact_master_problem(
                data,
                graph,
                artifacts,
                master_build_options
            );
            if (args.gurobi_threads >= 0) {
                problem.model->set(GRB_IntParam_Threads, args.gurobi_threads);
            }
            problem.model->set(GRB_DoubleParam_TimeLimit, args.solver_time_limit);

            output_file << "[main] Compact master model built successfully.\n";
            if (args.solve_model == 1) {
                spdp::solve_compact_master_problem(problem);
                print_solution_summary(output_file, problem);
                const spdp::RecoveredSolution recovered_solution =
                    spdp::recover_incumbent_solution(data, graph, artifacts, problem);
                spdp::write_recovered_solution(solution_file, recovered_solution);
            } else {
                output_file << "[main] Solve skipped by CLI option.\n";
                solution_file << "!! Print 2\n";
                solution_file << "Solution:\n";
                solution_file << "  Solve skipped by CLI option\n";
                solution_file << "Solution done \n";
            }
        }

        std::cout << "Log written to " << log_output_path << '\n';
        std::cout << "Solution written to " << solution_output_path << '\n';
        return 0;
    } catch (const GRBException& ex) {
        std::cerr << "Gurobi error (" << ex.getErrorCode() << "): " << ex.getMessage() << '\n';
        return 1;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        print_usage(argv[0]);
        return 1;
    }
}
