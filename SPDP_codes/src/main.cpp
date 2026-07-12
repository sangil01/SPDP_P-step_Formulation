#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
    std::string enumeration_sos1_mode = "default";
    std::string bnp_tree_mode = "root-only";
    std::string node_cg_phase_one_mode = "exact-cg";
    std::string node_cg_phase_two_pricing_mode = "exact-pricing";
    std::string node_cg_phase1_lp_method = "primal";
    std::string node_cg_phase2_lp_method = "primal";
    std::string bnp_branching_rule = "closest-to-half";
    std::string bnp_node_selection_rule = "best-bound";
    double bnp_theta_integrality_tolerance = 1e-6;
    double bnp_gap_tolerance = 1e-6;
    double bnp_initial_upper_bound = -1.0;
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
    int add_vi_36_combined = 0;
    int vi_36_subset_max_size = 1;
    int add_vi_request_block_sec = 0;
    int vi_request_block_sec_max_size = 2;
    int add_vi_44 = 0;
    int cg_max_iterations_per_phase = 1000;
    int exact_pricing_max_columns_per_start = 16;
    int exact_pricing_max_total_columns_per_round = 1024;
    double cg_reduced_cost_tolerance = -1e-6;
    int phase_one_heuristic_cg_max_attempts = 0;
    int phase_one_heuristic_cg_max_incumbents = 20;
    int phase_one_heuristic_cg_top_l = 3;
    double phase_one_heuristic_cg_weight_cost = 1.0;
    double phase_one_heuristic_cg_weight_time = 0.1;
    double phase_one_heuristic_cg_weight_saving = 0.5;
    int phase_one_heuristic_cg_random_seed = 1;
    int phase_two_heuristic_max_starts = 128;
    double phase_two_heuristic_start_ratio = 0.25;
    int phase_two_heuristic_ladder_levels = 3;
    int phase_two_heuristic_max_columns_per_start = 1;
    int phase_two_heuristic_max_columns_total = 256;
    double phase_two_heuristic_search_column_ratio = 1.0;
    std::string phase_two_heuristic_start_score_mode = "one-step-min";
    std::string phase_two_heuristic_engine = "shallow-search";
    int phase_two_labeling_top_k_next = 0;
    int phase_two_shallow_k1 = 8;
    int phase_two_shallow_k2 = 4;
    int phase_two_column_pool_enable = 0;
    int phase_two_column_pool_max_size = 5000;
    int phase_two_column_pool_max_reprice = 256;
    std::string full_enumeration_rc_update_mode = "sequential";
    std::string full_enumeration_parallel_stage1_backend = "custom";
    std::string full_enumeration_parallel_stage2_backend = "custom";
    int full_enumeration_rc_update_threads = 0;
    int full_enumeration_rc_detail_log = 1;
};

void print_usage(const char* executable) {
    std::cerr << "Usage: " << executable
              << " [instance] [--p N] [--solver-time-limit T] [--gurobi-threads N]"
              << " [--solver-mode enumeration|branch-and-price]"
              << " [--enumeration-sos1-mode default|sos1-auto|sos1-native]"
              << " [--bnp-tree-mode root-only|full-tree]"
              << " [--node-cg-phase1-mode exact-cg|heuristic-cg|heuristic-cg-3-step]"
              << " [--node-cg-phase2-pricing-mode exact-pricing|heuristic-pricing-then-exact|full-enumeration]"
              << " [--node-cg-phase1-lp-method automatic|primal|dual|barrier|concurrent]"
              << " [--node-cg-phase2-lp-method automatic|primal|dual|barrier|concurrent]"
              << " [--bnp-branching-rule closest-to-half]"
              << " [--bnp-node-selection-rule best-bound|dfs]"
              << " [--bnp-theta-integrality-tolerance T]"
              << " [--bnp-gap-tolerance T]"
              << " [--bnp-initial-upper-bound UB]"
              << " [--vi-formulation theta|x] [--dump-psteps N] [--validate-psteps 0|1] [--solve 0|1]"
              << " [--prune-infeasible-edges 0|1] [--prune-dominated-edges 0|1]"
              << " [--prune-symmetry-40 0|1] [--prune-symmetry-41 0|1]"
              << " [--prune-pickup-symmetry-43 0|1]"
              << " [--prune-delivery-symmetry-43 0|1]"
              << " [--add-vi-35 0|1]"
              << " [--add-vi-36-combined 0|1] [--vi-36-subset-max-size N]"
              << " [--add-vi-request-block-sec 0|1]"
              << " [--vi-request-block-sec-max-size N]"
              << " [--add-vi-44 0|1]"
              << " [--cg-max-iterations-per-phase N]"
              << " [--exact-pricing-max-columns-per-start N]"
              << " [--exact-pricing-max-total-columns-per-round N]"
              << " [--cg-reduced-cost-tolerance T]"
              << " [--phase1-heuristic-cg-max-attempts N]"
              << " [--phase1-heuristic-cg-max-incumbents N]"
              << " [--phase1-heuristic-cg-top-l N]"
              << " [--phase1-heuristic-cg-weight-cost W]"
              << " [--phase1-heuristic-cg-weight-time W]"
              << " [--phase1-heuristic-cg-weight-saving W]"
              << " [--phase1-heuristic-cg-random-seed N]"
              << " [--phase2-heuristic-max-starts N]"
              << " [--phase2-heuristic-start-ratio R]"
              << " [--phase2-heuristic-ladder-levels L]"
              << " [--phase2-heuristic-max-columns-per-start N]"
              << " [--phase2-heuristic-max-columns-total N]"
              << " [--phase2-heuristic-search-column-ratio R]"
              << " [--phase2-heuristic-start-score-mode one-step-min]"
              << " [--phase2-heuristic-engine labeling|shallow-search]"
              << " [--phase2-labeling-top-k-next N]"
              << " [--phase2-shallow-k1 N] [--phase2-shallow-k2 N]"
              << " [--phase2-column-pool-enable 0|1]"
              << " [--phase2-column-pool-max-size N]"
              << " [--phase2-column-pool-max-reprice N]"
              << " [--full-enumeration-rc-update-mode sequential|parallel]"
              << " [--full-enumeration-parallel-stage1-backend custom|onemkl]"
              << " [--full-enumeration-parallel-stage2-backend custom|onemkl]"
              << " [--full-enumeration-rc-update-threads N]"
              << " [--full-enumeration-rc-detail-log 0|1]\n";
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

int parse_binary_flag(const std::string& value, const std::string& field_name) {
    const int parsed = parse_int(value, field_name);
    if (parsed != 0 && parsed != 1) {
        throw std::runtime_error(field_name + " must be 0 or 1.");
    }
    return parsed;
}

std::string parse_vi_formulation(const std::string& value) {
    if (value == "theta" || value == "x") {
        return value;
    }
    throw std::runtime_error("Invalid value for --vi-formulation: " + value + " (expected theta or x)");
}

std::string parse_solver_mode(const std::string& value) {
    if (value == "enumeration" || value == "branch-and-price") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for --solver-mode: " + value +
        " (expected enumeration or branch-and-price)"
    );
}

std::string parse_enumeration_sos1_mode(const std::string& value) {
    if (value == "default" || value == "sos1-auto" || value == "sos1-native") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for --enumeration-sos1-mode: " + value +
        " (expected default, sos1-auto, or sos1-native)"
    );
}

std::string parse_bnp_tree_mode(const std::string& value) {
    if (value == "root-only" || value == "full-tree") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for --bnp-tree-mode: " + value +
        " (expected root-only or full-tree)"
    );
}

std::string parse_node_cg_phase_one_mode(const std::string& value) {
    if (value == "exact-cg" || value == "heuristic-cg" ||
        value == "heuristic-cg-3-step") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for node CG phase-one mode: " + value +
        " (expected exact-cg, heuristic-cg, or heuristic-cg-3-step)"
    );
}

std::string parse_node_cg_phase_two_pricing_mode(const std::string& value) {
    if (value == "exact-pricing" || value == "heuristic-pricing-then-exact" ||
        value == "full-enumeration") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for node CG phase-two pricing mode: " + value +
        " (expected exact-pricing, heuristic-pricing-then-exact, or full-enumeration)"
    );
}

std::string parse_phase_two_heuristic_start_score_mode(const std::string& value) {
    if (value == "one-step-min") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for --phase2-heuristic-start-score-mode: " + value +
        " (expected one-step-min)"
    );
}

std::string parse_phase_two_heuristic_engine(const std::string& value) {
    if (value == "labeling" || value == "shallow-search") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for --phase2-heuristic-engine: " + value +
        " (expected labeling or shallow-search)"
    );
}

std::string parse_full_enumeration_rc_update_mode(const std::string& value) {
    if (value == "sequential" || value == "parallel") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for --full-enumeration-rc-update-mode: " + value +
        " (expected sequential or parallel)"
    );
}

std::string parse_full_enumeration_parallel_backend(const std::string& value) {
    if (value == "custom" || value == "onemkl") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for full-enumeration parallel backend: " + value +
        " (expected custom or onemkl)"
    );
}

std::string parse_lp_method(const std::string& value) {
    if (value == "automatic" || value == "primal" || value == "dual" ||
        value == "barrier" || value == "concurrent") {
        return value;
    }
    throw std::runtime_error(
        "Invalid LP method: " + value +
        " (expected automatic, primal, dual, barrier, or concurrent)"
    );
}

std::string parse_bnp_branching_rule(const std::string& value) {
    if (value == "closest-to-half") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for --bnp-branching-rule: " + value +
        " (expected closest-to-half)"
    );
}

std::string parse_bnp_node_selection_rule(const std::string& value) {
    if (value == "best-bound" || value == "dfs") {
        return value;
    }
    throw std::runtime_error(
        "Invalid value for --bnp-node-selection-rule: " + value +
        " (expected best-bound or dfs)"
    );
}

spdp::GurobiLPMethod to_lp_method(const std::string& value) {
    if (value == "automatic") {
        return spdp::GurobiLPMethod::Automatic;
    }
    if (value == "primal") {
        return spdp::GurobiLPMethod::Primal;
    }
    if (value == "dual") {
        return spdp::GurobiLPMethod::Dual;
    }
    if (value == "barrier") {
        return spdp::GurobiLPMethod::Barrier;
    }
    if (value == "concurrent") {
        return spdp::GurobiLPMethod::Concurrent;
    }
    throw std::runtime_error("Unsupported LP method: " + value);
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

spdp::EnumerationSOS1Mode to_enumeration_sos1_mode(const std::string& value) {
    if (value == "default") {
        return spdp::EnumerationSOS1Mode::Default;
    }
    if (value == "sos1-auto") {
        return spdp::EnumerationSOS1Mode::Auto;
    }
    if (value == "sos1-native") {
        return spdp::EnumerationSOS1Mode::Native;
    }
    throw std::runtime_error("Unsupported enumeration SOS1 mode: " + value);
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

        if (arg == "--enumeration-sos1-mode") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--enumeration-sos1-mode requires a value.");
            }
            options.enumeration_sos1_mode =
                parse_enumeration_sos1_mode(argv[++idx]);
            continue;
        }

        if (arg == "--bnp-tree-mode") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--bnp-tree-mode requires a value.");
            }
            options.bnp_tree_mode = parse_bnp_tree_mode(argv[++idx]);
            continue;
        }

        if (arg == "--node-cg-phase1-mode") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.node_cg_phase_one_mode = parse_node_cg_phase_one_mode(argv[++idx]);
            continue;
        }

        if (arg == "--node-cg-phase2-pricing-mode") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--node-cg-phase2-pricing-mode requires a value.");
            }
            options.node_cg_phase_two_pricing_mode =
                parse_node_cg_phase_two_pricing_mode(argv[++idx]);
            continue;
        }

        if (arg == "--node-cg-phase1-lp-method") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.node_cg_phase1_lp_method = parse_lp_method(argv[++idx]);
            continue;
        }

        if (arg == "--node-cg-phase2-lp-method") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.node_cg_phase2_lp_method = parse_lp_method(argv[++idx]);
            continue;
        }

        if (arg == "--bnp-branching-rule") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.bnp_branching_rule = parse_bnp_branching_rule(argv[++idx]);
            continue;
        }

        if (arg == "--bnp-node-selection-rule") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.bnp_node_selection_rule =
                parse_bnp_node_selection_rule(argv[++idx]);
            continue;
        }

        if (arg == "--bnp-theta-integrality-tolerance") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.bnp_theta_integrality_tolerance =
                parse_double(argv[++idx], arg);
            if (options.bnp_theta_integrality_tolerance <= 0.0) {
                throw std::runtime_error(arg + " must be positive.");
            }
            continue;
        }

        if (arg == "--bnp-gap-tolerance") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.bnp_gap_tolerance = parse_double(argv[++idx], arg);
            if (options.bnp_gap_tolerance < 0.0) {
                throw std::runtime_error(arg + " must be nonnegative.");
            }
            continue;
        }

        if (arg == "--bnp-initial-upper-bound") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.bnp_initial_upper_bound = parse_double(argv[++idx], arg);
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

        if (arg == "--exact-pricing-max-columns-per-start") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.exact_pricing_max_columns_per_start =
                parse_int(argv[++idx], arg);
            if (options.exact_pricing_max_columns_per_start <= 0) {
                throw std::runtime_error(arg + " must be positive.");
            }
            continue;
        }

        if (arg == "--exact-pricing-max-total-columns-per-round") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.exact_pricing_max_total_columns_per_round =
                parse_int(argv[++idx], arg);
            if (options.exact_pricing_max_total_columns_per_round <= 0) {
                throw std::runtime_error(arg + " must be positive.");
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

        if (arg == "--phase1-heuristic-cg-max-attempts") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.phase_one_heuristic_cg_max_attempts = parse_int(argv[++idx], arg);
            if (options.phase_one_heuristic_cg_max_attempts < 0) {
                throw std::runtime_error(arg + " must be nonnegative.");
            }
            continue;
        }

        if (arg == "--phase1-heuristic-cg-max-incumbents") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.phase_one_heuristic_cg_max_incumbents = parse_int(argv[++idx], arg);
            if (options.phase_one_heuristic_cg_max_incumbents <= 0) {
                throw std::runtime_error(arg + " must be positive.");
            }
            continue;
        }

        if (arg == "--phase1-heuristic-cg-top-l") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.phase_one_heuristic_cg_top_l = parse_int(argv[++idx], arg);
            if (options.phase_one_heuristic_cg_top_l <= 0) {
                throw std::runtime_error(arg + " must be positive.");
            }
            continue;
        }

        if (arg == "--phase1-heuristic-cg-weight-cost") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.phase_one_heuristic_cg_weight_cost = parse_double(argv[++idx], arg);
            continue;
        }

        if (arg == "--phase1-heuristic-cg-weight-time") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.phase_one_heuristic_cg_weight_time = parse_double(argv[++idx], arg);
            continue;
        }

        if (arg == "--phase1-heuristic-cg-weight-saving") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.phase_one_heuristic_cg_weight_saving = parse_double(argv[++idx], arg);
            continue;
        }

        if (arg == "--phase1-heuristic-cg-random-seed") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value.");
            }
            options.phase_one_heuristic_cg_random_seed = parse_int(argv[++idx], arg);
            if (options.phase_one_heuristic_cg_random_seed < 0) {
                throw std::runtime_error(arg + " must be nonnegative.");
            }
            continue;
        }

        if (arg == "--phase2-heuristic-max-starts") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-heuristic-max-starts requires a value.");
            }
            options.phase_two_heuristic_max_starts =
                parse_int(argv[++idx], "--phase2-heuristic-max-starts");
            if (options.phase_two_heuristic_max_starts <= 0) {
                throw std::runtime_error("--phase2-heuristic-max-starts must be positive.");
            }
            continue;
        }

        if (arg == "--phase2-heuristic-start-ratio") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-heuristic-start-ratio requires a value.");
            }
            options.phase_two_heuristic_start_ratio =
                parse_double(argv[++idx], "--phase2-heuristic-start-ratio");
            if (options.phase_two_heuristic_start_ratio <= 0.0 ||
                options.phase_two_heuristic_start_ratio > 1.0) {
                throw std::runtime_error(
                    "--phase2-heuristic-start-ratio must lie in (0, 1]."
                );
            }
            continue;
        }

        if (arg == "--phase2-heuristic-ladder-levels") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-heuristic-ladder-levels requires a value.");
            }
            options.phase_two_heuristic_ladder_levels =
                parse_int(argv[++idx], "--phase2-heuristic-ladder-levels");
            if (options.phase_two_heuristic_ladder_levels < 0) {
                throw std::runtime_error(
                    "--phase2-heuristic-ladder-levels must be nonnegative."
                );
            }
            continue;
        }

        if (arg == "--phase2-heuristic-max-columns-per-start") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(
                    "--phase2-heuristic-max-columns-per-start requires a value."
                );
            }
            options.phase_two_heuristic_max_columns_per_start =
                parse_int(argv[++idx], "--phase2-heuristic-max-columns-per-start");
            if (options.phase_two_heuristic_max_columns_per_start <= 0) {
                throw std::runtime_error(
                    "--phase2-heuristic-max-columns-per-start must be positive."
                );
            }
            continue;
        }

        if (arg == "--phase2-heuristic-max-columns-total") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-heuristic-max-columns-total requires a value.");
            }
            options.phase_two_heuristic_max_columns_total =
                parse_int(argv[++idx], "--phase2-heuristic-max-columns-total");
            if (options.phase_two_heuristic_max_columns_total <= 0) {
                throw std::runtime_error(
                    "--phase2-heuristic-max-columns-total must be positive."
                );
            }
            continue;
        }

        if (arg == "--phase2-heuristic-search-column-ratio") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-heuristic-search-column-ratio requires a value.");
            }
            options.phase_two_heuristic_search_column_ratio =
                parse_double(argv[++idx], "--phase2-heuristic-search-column-ratio");
            if (options.phase_two_heuristic_search_column_ratio < 1.0) {
                throw std::runtime_error(
                    "--phase2-heuristic-search-column-ratio must be at least 1.0."
                );
            }
            continue;
        }

        if (arg == "--phase2-heuristic-start-score-mode") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(
                    "--phase2-heuristic-start-score-mode requires a value."
                );
            }
            options.phase_two_heuristic_start_score_mode =
                parse_phase_two_heuristic_start_score_mode(argv[++idx]);
            continue;
        }

        if (arg == "--phase2-heuristic-engine") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-heuristic-engine requires a value.");
            }
            options.phase_two_heuristic_engine =
                parse_phase_two_heuristic_engine(argv[++idx]);
            continue;
        }

        if (arg == "--phase2-labeling-top-k-next") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-labeling-top-k-next requires a value.");
            }
            options.phase_two_labeling_top_k_next =
                parse_int(argv[++idx], "--phase2-labeling-top-k-next");
            if (options.phase_two_labeling_top_k_next < 0) {
                throw std::runtime_error(
                    "--phase2-labeling-top-k-next must be nonnegative."
                );
            }
            continue;
        }

        if (arg == "--phase2-shallow-k1") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-shallow-k1 requires a value.");
            }
            options.phase_two_shallow_k1 =
                parse_int(argv[++idx], "--phase2-shallow-k1");
            if (options.phase_two_shallow_k1 <= 0) {
                throw std::runtime_error("--phase2-shallow-k1 must be positive.");
            }
            continue;
        }

        if (arg == "--phase2-shallow-k2") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-shallow-k2 requires a value.");
            }
            options.phase_two_shallow_k2 =
                parse_int(argv[++idx], "--phase2-shallow-k2");
            if (options.phase_two_shallow_k2 <= 0) {
                throw std::runtime_error("--phase2-shallow-k2 must be positive.");
            }
            continue;
        }

        if (arg == "--phase2-column-pool-enable") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-column-pool-enable requires a value.");
            }
            options.phase_two_column_pool_enable =
                parse_int(argv[++idx], "--phase2-column-pool-enable");
            if (options.phase_two_column_pool_enable != 0 &&
                options.phase_two_column_pool_enable != 1) {
                throw std::runtime_error("--phase2-column-pool-enable must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--phase2-column-pool-max-size") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-column-pool-max-size requires a value.");
            }
            options.phase_two_column_pool_max_size =
                parse_int(argv[++idx], "--phase2-column-pool-max-size");
            if (options.phase_two_column_pool_max_size <= 0) {
                throw std::runtime_error("--phase2-column-pool-max-size must be positive.");
            }
            continue;
        }

        if (arg == "--phase2-column-pool-max-reprice") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--phase2-column-pool-max-reprice requires a value.");
            }
            options.phase_two_column_pool_max_reprice =
                parse_int(argv[++idx], "--phase2-column-pool-max-reprice");
            if (options.phase_two_column_pool_max_reprice <= 0) {
                throw std::runtime_error(
                    "--phase2-column-pool-max-reprice must be positive."
                );
            }
            continue;
        }

        if (arg == "--full-enumeration-rc-update-mode") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(
                    "--full-enumeration-rc-update-mode requires a value."
                );
            }
            options.full_enumeration_rc_update_mode =
                parse_full_enumeration_rc_update_mode(argv[++idx]);
            continue;
        }

        if (arg == "--full-enumeration-parallel-stage1-backend") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(
                    "--full-enumeration-parallel-stage1-backend requires a value."
                );
            }
            options.full_enumeration_parallel_stage1_backend =
                parse_full_enumeration_parallel_backend(argv[++idx]);
            continue;
        }

        if (arg == "--full-enumeration-parallel-stage2-backend") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(
                    "--full-enumeration-parallel-stage2-backend requires a value."
                );
            }
            options.full_enumeration_parallel_stage2_backend =
                parse_full_enumeration_parallel_backend(argv[++idx]);
            continue;
        }

        if (arg == "--full-enumeration-rc-update-threads") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(
                    "--full-enumeration-rc-update-threads requires a value."
                );
            }
            options.full_enumeration_rc_update_threads =
                parse_int(argv[++idx], "--full-enumeration-rc-update-threads");
            if (options.full_enumeration_rc_update_threads < 0) {
                throw std::runtime_error(
                    "--full-enumeration-rc-update-threads must be nonnegative."
                );
            }
            continue;
        }

        if (arg == "--full-enumeration-rc-detail-log") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(
                    "--full-enumeration-rc-detail-log requires a value."
                );
            }
            options.full_enumeration_rc_detail_log =
                parse_int(argv[++idx], "--full-enumeration-rc-detail-log");
            if (options.full_enumeration_rc_detail_log != 0 &&
                options.full_enumeration_rc_detail_log != 1) {
                throw std::runtime_error(
                    "--full-enumeration-rc-detail-log must be 0 or 1."
                );
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

        if (arg == "--add-vi-36-combined") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--add-vi-36-combined requires a value.");
            }
            options.add_vi_36_combined =
                parse_int(argv[++idx], "--add-vi-36-combined");
            if (options.add_vi_36_combined != 0 && options.add_vi_36_combined != 1) {
                throw std::runtime_error("--add-vi-36-combined must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--vi-36-subset-max-size") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--vi-36-subset-max-size requires a value.");
            }
            options.vi_36_subset_max_size =
                parse_int(argv[++idx], "--vi-36-subset-max-size");
            if (options.vi_36_subset_max_size < 1) {
                throw std::runtime_error("--vi-36-subset-max-size must be at least 1.");
            }
            continue;
        }

        if (arg == "--add-vi-request-block-sec") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--add-vi-request-block-sec requires a value.");
            }
            options.add_vi_request_block_sec =
                parse_int(argv[++idx], "--add-vi-request-block-sec");
            if (options.add_vi_request_block_sec != 0 &&
                options.add_vi_request_block_sec != 1) {
                throw std::runtime_error("--add-vi-request-block-sec must be 0 or 1.");
            }
            continue;
        }

        if (arg == "--vi-request-block-sec-max-size") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(
                    "--vi-request-block-sec-max-size requires a value."
                );
            }
            options.vi_request_block_sec_max_size =
                parse_int(argv[++idx], "--vi-request-block-sec-max-size");
            if (options.vi_request_block_sec_max_size < 2) {
                throw std::runtime_error(
                    "--vi-request-block-sec-max-size must be at least 2."
                );
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
    out << "[main] Enumeration SOS1 mode: " << args.enumeration_sos1_mode << '\n';
    out << "[main] BnP tree mode: " << args.bnp_tree_mode << '\n';
    out << "[main] Node CG phase-1 mode: " << args.node_cg_phase_one_mode << '\n';
    out << "[main] Node CG phase-2 pricing mode: "
        << args.node_cg_phase_two_pricing_mode << '\n';
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
    out << "[main] add_vi_36_combined: " << args.add_vi_36_combined << '\n';
    out << "[main] vi_36_subset_max_size: " << args.vi_36_subset_max_size << '\n';
    out << "[main] add_vi_request_block_sec: " << args.add_vi_request_block_sec << '\n';
    out << "[main] vi_request_block_sec_max_size: "
        << args.vi_request_block_sec_max_size << '\n';
    out << "[main] add_vi_44: " << args.add_vi_44 << '\n';
    out << "[main] cg_max_iterations_per_phase: " << args.cg_max_iterations_per_phase << '\n';
    out << "[main] exact_pricing_max_columns_per_start: "
        << args.exact_pricing_max_columns_per_start << '\n';
    out << "[main] exact_pricing_max_total_columns_per_round: "
        << args.exact_pricing_max_total_columns_per_round << '\n';
    out << "[main] phase1_heuristic_cg_max_attempts: "
        << args.phase_one_heuristic_cg_max_attempts << '\n';
    out << "[main] phase1_heuristic_cg_max_incumbents: "
        << args.phase_one_heuristic_cg_max_incumbents << '\n';
    out << "[main] phase1_heuristic_cg_top_l: "
        << args.phase_one_heuristic_cg_top_l << '\n';
    out << "[main] phase1_heuristic_cg_weight_cost: "
        << format_double(args.phase_one_heuristic_cg_weight_cost) << '\n';
    out << "[main] phase1_heuristic_cg_weight_time: "
        << format_double(args.phase_one_heuristic_cg_weight_time) << '\n';
    out << "[main] phase1_heuristic_cg_weight_saving: "
        << format_double(args.phase_one_heuristic_cg_weight_saving) << '\n';
    out << "[main] phase1_heuristic_cg_random_seed: "
        << args.phase_one_heuristic_cg_random_seed << '\n';
    out << "[main] phase2_heuristic_max_starts: " << args.phase_two_heuristic_max_starts << '\n';
    out << "[main] phase2_heuristic_start_ratio: "
        << format_double(args.phase_two_heuristic_start_ratio) << '\n';
    out << "[main] phase2_heuristic_ladder_levels: "
        << args.phase_two_heuristic_ladder_levels << '\n';
    out << "[main] phase2_heuristic_max_columns_per_start: "
        << args.phase_two_heuristic_max_columns_per_start << '\n';
    out << "[main] phase2_heuristic_max_columns_total: "
        << args.phase_two_heuristic_max_columns_total << '\n';
    out << "[main] phase2_heuristic_search_column_ratio: "
        << format_double(args.phase_two_heuristic_search_column_ratio) << '\n';
    out << "[main] phase2_heuristic_start_score_mode: "
        << args.phase_two_heuristic_start_score_mode << '\n';
    out << "[main] phase2_heuristic_engine: "
        << args.phase_two_heuristic_engine << '\n';
    out << "[main] phase2_labeling_top_k_next: "
        << args.phase_two_labeling_top_k_next << '\n';
    out << "[main] phase2_shallow_k1: " << args.phase_two_shallow_k1 << '\n';
    out << "[main] phase2_shallow_k2: " << args.phase_two_shallow_k2 << '\n';
    out << "[main] phase2_column_pool_enable: "
        << args.phase_two_column_pool_enable << '\n';
    out << "[main] phase2_column_pool_max_size: "
        << args.phase_two_column_pool_max_size << '\n';
    out << "[main] phase2_column_pool_max_reprice: "
        << args.phase_two_column_pool_max_reprice << '\n';
    out << "[main] full_enumeration_rc_update_mode: "
        << args.full_enumeration_rc_update_mode << '\n';
    out << "[main] full_enumeration_parallel_stage1_backend: "
        << args.full_enumeration_parallel_stage1_backend << '\n';
    out << "[main] full_enumeration_parallel_stage2_backend: "
        << args.full_enumeration_parallel_stage2_backend << '\n';
    out << "[main] full_enumeration_rc_update_threads: "
        << args.full_enumeration_rc_update_threads << '\n';
    out << "[main] full_enumeration_rc_detail_log: "
        << args.full_enumeration_rc_detail_log << '\n';
    out << "[main] node_cg_phase1_lp_method: "
        << args.node_cg_phase1_lp_method << '\n';
    out << "[main] node_cg_phase2_lp_method: "
        << args.node_cg_phase2_lp_method << '\n';
    out << "[main] bnp_branching_rule: "
        << args.bnp_branching_rule << '\n';
    out << "[main] bnp_node_selection_rule: "
        << args.bnp_node_selection_rule << '\n';
    out << "[main] bnp_theta_integrality_tolerance: "
        << format_double(args.bnp_theta_integrality_tolerance) << '\n';
    out << "[main] bnp_gap_tolerance: "
        << format_double(args.bnp_gap_tolerance) << '\n';
    out << "[main] bnp_initial_upper_bound: "
        << format_double(args.bnp_initial_upper_bound) << '\n';
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

        if (args.solver_mode == "branch-and-price") {
            if (args.vi_formulation != "theta") {
                throw std::runtime_error(
                    "branch-and-price mode currently supports only --vi-formulation theta."
                );
            }

            if (args.solve_model == 1) {
                spdp::NodeCGOptions cg_options;
                cg_options.p = args.p;
                cg_options.time_limit = data.time_limit;
                cg_options.solver_time_limit = args.solver_time_limit;
                cg_options.gurobi_threads = args.gurobi_threads;
                cg_options.phase_one_lp_method =
                    to_lp_method(args.node_cg_phase1_lp_method);
                cg_options.phase_two_lp_method =
                    to_lp_method(args.node_cg_phase2_lp_method);
                cg_options.max_iterations_per_phase =
                    static_cast<std::size_t>(args.cg_max_iterations_per_phase);
                cg_options.exact_pricing_max_columns_per_start =
                    static_cast<std::size_t>(args.exact_pricing_max_columns_per_start);
                cg_options.exact_pricing_max_total_columns_per_round =
                    static_cast<std::size_t>(args.exact_pricing_max_total_columns_per_round);
                cg_options.reduced_cost_tolerance = args.cg_reduced_cost_tolerance;
                cg_options.phase_one_heuristic_cg_max_attempts =
                    static_cast<std::size_t>(args.phase_one_heuristic_cg_max_attempts);
                cg_options.phase_one_heuristic_cg_max_incumbents =
                    static_cast<std::size_t>(args.phase_one_heuristic_cg_max_incumbents);
                cg_options.phase_one_heuristic_cg_top_l =
                    static_cast<std::size_t>(args.phase_one_heuristic_cg_top_l);
                cg_options.phase_one_heuristic_cg_weight_cost =
                    args.phase_one_heuristic_cg_weight_cost;
                cg_options.phase_one_heuristic_cg_weight_time =
                    args.phase_one_heuristic_cg_weight_time;
                cg_options.phase_one_heuristic_cg_weight_saving =
                    args.phase_one_heuristic_cg_weight_saving;
                cg_options.phase_one_heuristic_cg_random_seed =
                    static_cast<unsigned int>(args.phase_one_heuristic_cg_random_seed);
                cg_options.prune_pickup_symmetry_43 = args.prune_pickup_symmetry_43 == 1;
                cg_options.prune_delivery_symmetry_43 = args.prune_delivery_symmetry_43 == 1;
                cg_options.add_root_vi_35 = args.add_vi_35 == 1;
                cg_options.add_root_vi_36_combined = args.add_vi_36_combined == 1;
                cg_options.root_vi_36_subset_max_size =
                    static_cast<std::size_t>(args.vi_36_subset_max_size);
                cg_options.add_root_vi_request_block_sec =
                    args.add_vi_request_block_sec == 1;
                cg_options.root_vi_request_block_sec_max_size =
                    static_cast<std::size_t>(args.vi_request_block_sec_max_size);
                cg_options.add_root_vi_44 = args.add_vi_44 == 1;
                cg_options.gurobi_log_path = gurobi_log_path.string();
                if (args.node_cg_phase_one_mode == "exact-cg") {
                    cg_options.phase_one_mode = spdp::NodeCGPhaseOneMode::ExactCG;
                } else if (args.node_cg_phase_one_mode == "heuristic-cg-3-step") {
                    cg_options.phase_one_mode = spdp::NodeCGPhaseOneMode::HeuristicCG3Step;
                } else {
                    cg_options.phase_one_mode = spdp::NodeCGPhaseOneMode::HeuristicCG;
                }
                if (args.node_cg_phase_two_pricing_mode == "exact-pricing") {
                    cg_options.phase_two_pricing_mode =
                        spdp::NodeCGPhaseTwoPricingMode::ExactPricing;
                } else if (args.node_cg_phase_two_pricing_mode == "full-enumeration") {
                    cg_options.phase_two_pricing_mode =
                        spdp::NodeCGPhaseTwoPricingMode::FullEnumeration;
                } else {
                    cg_options.phase_two_pricing_mode =
                        spdp::NodeCGPhaseTwoPricingMode::HeuristicPricingThenExact;
                }
                cg_options.phase_two_heuristic_max_starts =
                    static_cast<std::size_t>(args.phase_two_heuristic_max_starts);
                cg_options.phase_two_heuristic_start_ratio =
                    args.phase_two_heuristic_start_ratio;
                cg_options.phase_two_heuristic_ladder_levels =
                    static_cast<std::size_t>(args.phase_two_heuristic_ladder_levels);
                cg_options.phase_two_heuristic_max_columns_per_start =
                    static_cast<std::size_t>(args.phase_two_heuristic_max_columns_per_start);
                cg_options.phase_two_heuristic_max_total_columns =
                    static_cast<std::size_t>(args.phase_two_heuristic_max_columns_total);
                cg_options.phase_two_heuristic_search_column_ratio =
                    args.phase_two_heuristic_search_column_ratio;
                cg_options.phase_two_heuristic_start_score_mode =
                    spdp::HeuristicStartScoreMode::OneStepMin;
                cg_options.phase_two_heuristic_engine =
                    args.phase_two_heuristic_engine == "labeling"
                        ? spdp::NodeCGPhaseTwoHeuristicEngine::Labeling
                        : spdp::NodeCGPhaseTwoHeuristicEngine::ShallowSearch;
                cg_options.phase_two_labeling_top_k_next =
                    static_cast<std::size_t>(args.phase_two_labeling_top_k_next);
                cg_options.phase_two_shallow_k1 =
                    static_cast<std::size_t>(args.phase_two_shallow_k1);
                cg_options.phase_two_shallow_k2 =
                    static_cast<std::size_t>(args.phase_two_shallow_k2);
                cg_options.phase_two_column_pool_enabled =
                    args.phase_two_column_pool_enable == 1;
                cg_options.phase_two_column_pool_max_size =
                    static_cast<std::size_t>(args.phase_two_column_pool_max_size);
                cg_options.phase_two_column_pool_max_reprice =
                    static_cast<std::size_t>(args.phase_two_column_pool_max_reprice);
                cg_options.full_enumeration_rc_update_mode =
                    args.full_enumeration_rc_update_mode == "parallel"
                        ? spdp::FullEnumerationRCUpdateMode::Parallel
                        : spdp::FullEnumerationRCUpdateMode::Sequential;
                cg_options.full_enumeration_parallel_stage1_backend =
                    args.full_enumeration_parallel_stage1_backend == "onemkl"
                        ? spdp::FullEnumerationRCUpdateBackend::OneMKL
                        : spdp::FullEnumerationRCUpdateBackend::Custom;
                cg_options.full_enumeration_parallel_stage2_backend =
                    args.full_enumeration_parallel_stage2_backend == "onemkl"
                        ? spdp::FullEnumerationRCUpdateBackend::OneMKL
                        : spdp::FullEnumerationRCUpdateBackend::Custom;
                cg_options.full_enumeration_rc_update_threads =
                    static_cast<std::size_t>(args.full_enumeration_rc_update_threads);
                cg_options.full_enumeration_rc_detail_log =
                    args.full_enumeration_rc_detail_log == 1;
                if (args.bnp_tree_mode == "root-only") {
                    const spdp::NodeCGResult cg_result =
                        spdp::solve_node_column_generation(data, graph, cg_options, &output_file);
                    spdp::write_node_cg_summary(output_file, cg_result);
                    spdp::write_node_cg_solution(solution_file, cg_result, graph);
                } else {
                    spdp::BranchAndPriceOptions bnp_options;
                    bnp_options.node_cg_options = cg_options;
                    bnp_options.tree_mode = spdp::BnPTreeMode::FullTree;
                    bnp_options.branching_rule =
                        args.bnp_branching_rule == "closest-to-half"
                            ? spdp::BnPBranchingRule::ClosestToHalf
                            : spdp::BnPBranchingRule::ClosestToHalf;
                    bnp_options.node_selection_rule =
                        args.bnp_node_selection_rule == "dfs"
                            ? spdp::BnPNodeSelectionRule::DepthFirst
                            : spdp::BnPNodeSelectionRule::BestBound;
                    bnp_options.theta_integrality_tolerance =
                        args.bnp_theta_integrality_tolerance;
                    bnp_options.gap_tolerance = args.bnp_gap_tolerance;
                    bnp_options.instance_name =
                        std::filesystem::path(args.instance).stem().string();
                    bnp_options.initial_upper_bound = args.bnp_initial_upper_bound;
                    const spdp::BranchAndPriceResult bnp_result =
                        spdp::solve_branch_and_price(data, graph, bnp_options, &output_file);
                    spdp::write_branch_and_price_summary(output_file, bnp_result);
                    spdp::write_branch_and_price_solution(solution_file, bnp_result, graph);
                }
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
                args.add_vi_36_combined == 1,
                static_cast<std::size_t>(args.vi_36_subset_max_size),
                args.add_vi_request_block_sec == 1,
                static_cast<std::size_t>(args.vi_request_block_sec_max_size),
                args.add_vi_44 == 1,
                to_vi_formulation(args.vi_formulation),
                to_enumeration_sos1_mode(args.enumeration_sos1_mode),
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
            output_file << "[main] Parallel-edge SOS1 count: "
                        << problem.parallel_edge_sos1_count << '\n';
            output_file << "[main] Parallel-edge SOS1 max size: "
                        << problem.parallel_edge_sos1_max_size << '\n';
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
