#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <set>

#include "gurobi_c++.h"
#include "GenMultiGraph.h"
#include "ReadData.h"

namespace {

struct CliOptions {
    std::string instance = "A3.dat";
    int sample_edges = 5;
    int prune_infeasible_edges = 1;
    int prune_dominated_edges = 1;
};

void print_usage(const char* executable) {
    std::cerr << "Usage: " << executable
              << " [instance] [--sample-edges N] [--prune-infeasible-edges 0|1]"
              << " [--prune-dominated-edges 0|1]\n";
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

CliOptions parse_cli(int argc, char** argv) {
    CliOptions options;
    bool instance_set = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            std::exit(0);
        }

        if (arg == "--sample-edges") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--sample-edges requires a value.");
            }
            options.sample_edges = parse_int(argv[++i], "--sample-edges");
            continue;
        }

        if (arg == "--prune-dominated-edges") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--prune-dominated-edges requires a value.");
            }
            options.prune_dominated_edges = parse_int(argv[++i], "--prune-dominated-edges");
            if (options.prune_dominated_edges != 0 && options.prune_dominated_edges != 1) {
                throw std::runtime_error("--prune-dominated-edges must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--prune-infeasible-edges") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--prune-infeasible-edges requires a value.");
            }
            options.prune_infeasible_edges = parse_int(argv[++i], "--prune-infeasible-edges");
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

std::string format_sequence_pi(const std::vector<int>& sequence_pi) {
    if (sequence_pi.empty()) {
        return "()";
    }

    std::ostringstream out;
    out << "(";
    for (std::size_t i = 0; i < sequence_pi.size(); ++i) {
        if (i > 0) {
            out << ", ";
        }
        out << sequence_pi[i];
    }
    if (sequence_pi.size() == 1U) {
        out << ",";
    }
    out << ")";

    return out.str();
}

std::string format_two_decimals(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

std::filesystem::path project_root_path() {
    const std::filesystem::path source_path(__FILE__);
    return source_path.parent_path().parent_path().parent_path();
}

std::filesystem::path build_output_path(const std::string& instance) {
    const std::filesystem::path instance_path(instance);
    const std::string output_name = instance_path.stem().string() + ".txt";
    return project_root_path() / "SPDP_output" / output_name;
}

void print_data_info(
    std::ostream& out,
    const CliOptions& args,
    const spdp::SPDPData& data,
    const spdp::MultiDiGraph& graph
) {
    out << "[main] Loaded instance: " << args.instance << '\n';
    out << "[main] Locations: " << data.locations << '\n';
    out << "[main] Fixed vehicle cost: " << data.fixed_vehicle_cost << '\n';
    out << "[main] Time pick-up: " << data.time_pickup << '\n';
    out << "[main] Time empty: " << data.time_empty << '\n';
    out << "[main] Time delivery: " << data.time_delivery << '\n';
    out << "[main] Time limit: " << data.time_limit << '\n';
    out << "[main] Requests: " << data.requests.size() << '\n';
    out << "[main] Time matrix size: " << data.time.size() << "x"
        << (data.time.empty() ? 0 : data.time[0].size()) << '\n';
    out << "[main] Distance matrix size: " << data.distance.size() << "x"
        << (data.distance.empty() ? 0 : data.distance[0].size()) << '\n';
    out << "[main] Node count: " << graph.number_of_nodes() << '\n';
    out << "[main] Edge count: " << graph.number_of_edges() << '\n';
}

void print_selected_edge_info(
    std::ostream& out,
    const spdp::MultiDiGraph& graph,
    const std::vector<std::pair<int, int>>& target_pairs
) {
    out << "[main] Full edge list for selected node pairs:\n";

    for (const auto& pair : target_pairs) {
        const int u = pair.first;
        const int v = pair.second;

        const std::vector<spdp::EdgeRecord> pair_edges = graph.edges_between(u, v);
        out << "  Pair (" << u << ", " << v << ") -> " << pair_edges.size()
            << " edges\n";

        for (const spdp::EdgeRecord& edge : pair_edges) {
            out << "    key=" << edge.key
                << " | pi=" << format_sequence_pi(edge.data.sequence_pi)
                << " | time=" << format_two_decimals(edge.data.time)
                << " | cost=" << format_two_decimals(edge.data.cost)
                << " | start=" << spdp::state_to_str(edge.data.start_state)
                << " | end=" << spdp::state_to_str(edge.data.end_state)
                << '\n';
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const CliOptions args = parse_cli(argc, argv);
        const std::filesystem::path output_path = build_output_path(args.instance);
        std::filesystem::create_directories(output_path.parent_path());

        std::ofstream output_file(output_path);
        if (!output_file) {
            throw std::runtime_error("Failed to open output file: " + output_path.string());
        }

        const spdp::SPDPData data = spdp::read_spdp_data(args.instance);
        const spdp::MultiDiGraph graph = spdp::build_multigraph(
            data,
            args.prune_infeasible_edges == 1,
            args.prune_dominated_edges == 1,
            &output_file
        );

        print_data_info(output_file, args, data, graph);

        std::set<std::pair<int, int>> unique_pairs;
        for (const spdp::EdgeRecord& edge : graph.edges()) {
            unique_pairs.insert({edge.u, edge.v});
        }
        const std::vector<std::pair<int, int>> target_pairs(
            unique_pairs.begin(),
            unique_pairs.end()
        );

        print_selected_edge_info(output_file, graph, target_pairs);

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        print_usage(argv[0]);
        return 1;
    }
}
