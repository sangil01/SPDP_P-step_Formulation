#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "GenMultiGraph.h"
#include "ReadData.h"

namespace {

struct CliOptions {
    std::string instance = "A0.dat";
    int sample_edges = 5;
    int prune_dominated_edges = 0;
};

void print_usage(const char* executable) {
    std::cerr << "Usage: " << executable
              << " [instance] [--sample-edges N] [--prune-dominated-edges 0|1]\n";
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

}  // namespace

int main(int argc, char** argv) {
    try {
        const CliOptions args = parse_cli(argc, argv);

        const spdp::SPDPData data = spdp::read_spdp_data(args.instance);
        const spdp::MultiDiGraph graph = spdp::build_multigraph(
            data,
            args.prune_dominated_edges == 1
        );

        std::cout << "[main] Loaded instance: " << args.instance << '\n';
        std::cout << "[main] Requests: " << data.requests.size() << '\n';
        std::cout << "[main] Locations: " << data.locations << '\n';
        std::cout << "[main] Node count: " << graph.number_of_nodes() << '\n';
        std::cout << "[main] Edge count: " << graph.number_of_edges() << '\n';

        if (args.sample_edges > 0) {
            std::cout << "[main] Sample edges (up to " << args.sample_edges << "):\n";

            int shown = 0;
            for (const spdp::EdgeRecord& edge : graph.edges()) {
                std::cout << "  " << edge.u << " -> " << edge.v
                          << " | pi=" << format_sequence_pi(edge.data.sequence_pi)
                          << " | time=" << format_two_decimals(edge.data.time)
                          << " | cost=" << format_two_decimals(edge.data.cost)
                          << '\n';

                ++shown;
                if (shown >= args.sample_edges) {
                    break;
                }
            }
        }

        const std::vector<std::pair<int, int>> target_pairs = {
            {0, 1}, {1, 2}, {1, 3}, {3, 1}, {3, 4},
        };

        std::cout << "[main] Full edge list for selected node pairs:\n";

        for (const auto& pair : target_pairs) {
            const int u = pair.first;
            const int v = pair.second;

            const std::vector<spdp::EdgeRecord> pair_edges = graph.edges_between(u, v);
            std::cout << "  Pair (" << u << ", " << v << ") -> " << pair_edges.size()
                      << " edges\n";

            for (const spdp::EdgeRecord& edge : pair_edges) {
                std::cout << "    key=" << edge.key
                          << " | pi=" << format_sequence_pi(edge.data.sequence_pi)
                          << " | time=" << format_two_decimals(edge.data.time)
                          << " | cost=" << format_two_decimals(edge.data.cost)
                          << " | start=" << edge.data.start_state_str
                          << " | end=" << edge.data.end_state_str
                          << '\n';
            }
        }

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        print_usage(argv[0]);
        return 1;
    }
}
