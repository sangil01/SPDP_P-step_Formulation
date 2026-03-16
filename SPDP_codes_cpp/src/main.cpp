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
#include "PstepFormulation.h"
#include "ReadData.h"

namespace {

struct CliOptions {
    std::string instance = "A0.dat";
    int p = 2;
    double solver_time_limit = 3600.0;
    int dump_psteps = 10;
    int validate_psteps = 0;
    int solve_model = 1;
    int prune_infeasible_edges = 1;
    int prune_dominated_edges = 1;
};

void print_usage(const char* executable) {
    std::cerr << "Usage: " << executable
              << " [instance] [--p N] [--solver-time-limit T]"
              << " [--dump-psteps N] [--validate-psteps 0|1] [--solve 0|1]"
              << " [--prune-infeasible-edges 0|1] [--prune-dominated-edges 0|1]\n";
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

        if (arg == "--dump-psteps") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--dump-psteps requires a value.");
            }
            options.dump_psteps = parse_int(argv[++idx], "--dump-psteps");
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
    out << "[main] Solver time limit: " << format_double(args.solver_time_limit) << '\n';
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

        const spdp::MultiDiGraph graph = spdp::build_multigraph(
            data,
            args.prune_infeasible_edges == 1,
            args.prune_dominated_edges == 1,
            &output_file
        );

        print_instance_summary(output_file, args, data, graph);
        print_selected_edge_info(output_file, graph);

        const spdp::CompactPStepOptions options{
            args.p,
            data.time_limit,
            static_cast<std::size_t>(std::max(args.dump_psteps, 0)),
            args.validate_psteps == 1,
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

        spdp::CompactMasterProblem problem = spdp::build_compact_master_problem(
            graph,
            artifacts,
            build_gurobi_log_path(args.instance, args.p).string()
        );
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
