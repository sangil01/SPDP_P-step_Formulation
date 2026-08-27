#include "TestSupport.h"

#include "CpSatSpdpSolver.h"
#include "CpIncumbentAdapter.h"
#include "GenMultiGraph.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace {

spdp::SPDPData one_request_data() {
    spdp::SPDPData data;
    data.time_pickup = 2;
    data.time_empty = 3;
    data.time_delivery = 2;
    data.time_limit = 30;
    data.locations = 2;
    data.requests.push_back(spdp::Request{0, 1, 4});
    data.time = {{0, 4}, {4, 0}};
    data.distance = data.time;
    return data;
}

}  // namespace

SPDP_TEST(cp_data) {
    auto data = one_request_data();
    spdp::CpSatSolveOptions options;
    SPDP_CHECK(options.redundant.terminal_balance);
    SPDP_CHECK(options.redundant.full_skip_reservoir);
    SPDP_CHECK(!options.redundant.container_workload);
    SPDP_CHECK(options.redundant.aggregate_duration);
    SPDP_CHECK(!options.symmetry.first_pickup_vehicle_ordering);
    SPDP_CHECK(!options.symmetry.cor_43_identical_pickup_time_ordering);
    options.vehicle_count = 1;
    options.time_limit_seconds = 5;

    data.time[0][1] = 1.5;
    const auto invalid = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(invalid.outcome, spdp::CpSolveOutcome::ModelInvalid);

    data = one_request_data();
    const auto feasible = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(feasible.outcome, spdp::CpSolveOutcome::Feasible);
    SPDP_CHECK_EQ(feasible.routes.size(), 1U);
    SPDP_CHECK_EQ(feasible.routes.front().actions.size(), 3U);

    options.vehicle_count = 2;
    const auto infeasible = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(infeasible.outcome, spdp::CpSolveOutcome::ProvenInfeasible);
}

SPDP_TEST(cp_log_file) {
    const std::filesystem::path log_path =
        std::filesystem::temp_directory_path() / "spdp_cp_sat_progress.log";
    std::filesystem::remove(log_path);

    auto data = one_request_data();
    spdp::CpSatSolveOptions options;
    options.vehicle_count = 1;
    options.time_limit_seconds = 5;
    options.log_file_path = log_path.string();

    const auto solved = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(solved.outcome, spdp::CpSolveOutcome::Feasible);
    SPDP_CHECK(std::filesystem::exists(log_path));

    std::ifstream input(log_path);
    std::ostringstream contents;
    contents << input.rdbuf();
    SPDP_CHECK(contents.str().find("CP-SAT") != std::string::npos);
    SPDP_CHECK(contents.str().find("CpSolverResponse summary") !=
               std::string::npos);
    std::filesystem::remove(log_path);

    const std::filesystem::path invalid_log_path =
        std::filesystem::temp_directory_path() /
        "spdp_cp_sat_model_invalid.log";
    std::filesystem::remove(invalid_log_path);
    data.time[0][1] = 1.5;
    options.log_file_path = invalid_log_path.string();
    const auto invalid = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(invalid.outcome, spdp::CpSolveOutcome::ModelInvalid);
    SPDP_CHECK(std::filesystem::exists(invalid_log_path));
    std::ifstream invalid_input(invalid_log_path);
    std::ostringstream invalid_contents;
    invalid_contents << invalid_input.rdbuf();
    SPDP_CHECK(invalid_contents.str().find("MODEL_INVALID") !=
               std::string::npos);
    std::filesystem::remove(invalid_log_path);
}

SPDP_TEST(cp_symmetry) {
    auto data = one_request_data();
    data.time_limit = 100;
    data.locations = 7;
    data.requests = {
        spdp::Request{1, 5, 2},
        spdp::Request{1, 5, 2},
        spdp::Request{1, 6, 3},
    };
    data.time.assign(7, std::vector<double>(7, 1));
    for (int i = 0; i < 7; ++i) {
        data.time[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 0;
    }
    data.distance = data.time;

    spdp::CpSatSolveOptions options;
    options.vehicle_count = 1;
    options.time_limit_seconds = 5;
    const auto base = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(base.outcome, spdp::CpSolveOutcome::Feasible);
    SPDP_CHECK_EQ(base.build_stats.core_precedence_pruned_arcs, 3);
    SPDP_CHECK(base.build_stats.cor_40_pruned_arcs > 0);
    SPDP_CHECK(base.build_stats.cor_41_pruned_arcs > 0);
    SPDP_CHECK_EQ(base.build_stats.cor_43_constraint_count, 0);

    options.symmetry.cor_43_identical_pickup_time_ordering = true;
    const auto ordered = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(ordered.outcome, spdp::CpSolveOutcome::Feasible);
    SPDP_CHECK_EQ(ordered.build_stats.cor_43_constraint_count, 1);
    SPDP_CHECK_EQ(
        ordered.build_stats.created_action_arcs,
        base.build_stats.created_action_arcs
    );
}

SPDP_TEST(cp_core) {
    spdp::SPDPData data;
    data.time_pickup = 2;
    data.time_empty = 3;
    data.time_delivery = 2;
    data.time_limit = 100;
    data.locations = 4;
    data.requests = {
        spdp::Request{0, 3, 1},
        spdp::Request{1, 3, 2},
        spdp::Request{2, 3, 1},
    };
    data.time.assign(4, std::vector<double>(4, 1));
    for (int i = 0; i < 4; ++i) {
        data.time[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 0;
    }
    data.distance = data.time;

    spdp::CpSatSolveOptions options;
    options.vehicle_count = 1;
    options.time_limit_seconds = 5;
    const auto solved = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(solved.outcome, spdp::CpSolveOutcome::Feasible);
    SPDP_CHECK_EQ(solved.routes.size(), 1U);

    int load = 0;
    std::map<int, int> empty_inventory;
    std::vector<bool> picked(data.requests.size(), false);
    std::vector<bool> treated(data.requests.size(), false);
    const auto& route = solved.routes.front();
    SPDP_CHECK_EQ(route.actions.front().kind, spdp::CpActionKind::Pickup);
    SPDP_CHECK_EQ(route.actions.back().kind, spdp::CpActionKind::Delivery);
    SPDP_CHECK(route.completion_time <= data.time_limit);
    for (const auto& action : route.actions) {
        const auto& request =
            data.requests[static_cast<std::size_t>(action.request_index)];
        if (action.kind == spdp::CpActionKind::Pickup) {
            picked[static_cast<std::size_t>(action.request_index)] = true;
            ++load;
            SPDP_CHECK(load <= 2);
        } else if (action.kind == spdp::CpActionKind::Treatment) {
            SPDP_CHECK(picked[static_cast<std::size_t>(action.request_index)]);
            SPDP_CHECK(!treated[static_cast<std::size_t>(action.request_index)]);
            treated[static_cast<std::size_t>(action.request_index)] = true;
            ++empty_inventory[request.container_type];
        } else {
            SPDP_CHECK(empty_inventory[request.container_type] > 0);
            --empty_inventory[request.container_type];
            --load;
            SPDP_CHECK(load >= 0);
        }
    }
    SPDP_CHECK_EQ(load, 0);
    for (const auto& entry : empty_inventory) {
        SPDP_CHECK_EQ(entry.second, 0);
    }
}

SPDP_TEST(benchmark_a6) {
    const char* data_dir = std::getenv("SPDP_TEST_DATA_DIR");
    SPDP_CHECK(data_dir != nullptr);
    const auto data = spdp::read_spdp_data(
        std::string(data_dir) + "/RecDep_day_A6.dat"
    );
    const auto graph = spdp::build_multigraph(data);

    spdp::CpSatSolveOptions options;
    options.time_limit_seconds = 1200;
    options.workers = 0;
    options.vehicle_count = 1;
    const auto k1 = spdp::solve_fixed_k_cp_sat(data, options);
    std::cerr << "A6 K=1 status=" << k1.status_name
              << " raw=" << k1.raw_status
              << " error=" << k1.error_message
              << " wall=" << k1.wall_time_seconds << '\n';
    SPDP_CHECK_EQ(k1.outcome, spdp::CpSolveOutcome::ProvenInfeasible);

    options.vehicle_count = 2;
    const auto k2 = spdp::solve_fixed_k_cp_sat(data, options);
    std::cerr << "A6 K=2 status=" << k2.status_name
              << " raw=" << k2.raw_status
              << " error=" << k2.error_message
              << " wall=" << k2.wall_time_seconds << '\n';
    SPDP_CHECK_EQ(k2.outcome, spdp::CpSolveOutcome::Feasible);
    const auto mapped = spdp::map_cp_incumbent_to_multigraph(
        data, graph, k2.routes
    );
    SPDP_CHECK(mapped.success);
}

SPDP_TEST(benchmark_a6_k2_default) {
    const char* data_dir = std::getenv("SPDP_TEST_DATA_DIR");
    SPDP_CHECK(data_dir != nullptr);
    const auto data = spdp::read_spdp_data(
        std::string(data_dir) + "/RecDep_day_A6.dat"
    );
    const auto graph = spdp::build_multigraph(data);
    spdp::CpSatSolveOptions options;
    options.vehicle_count = 2;
    options.time_limit_seconds = 20;
    options.workers = 0;
    const auto cp = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK_EQ(cp.outcome, spdp::CpSolveOutcome::Feasible);
    const auto mapped = spdp::map_cp_incumbent_to_multigraph(
        data, graph, cp.routes
    );
    std::cerr << "A6 default K=2 adapter=" << mapped.error_message << '\n';
    SPDP_CHECK(mapped.success);
}

namespace {

void run_known_instance_case(
    const std::string& filename,
    int vehicle_count,
    spdp::CpSolveOutcome expected,
    bool require_adapter
) {
    const char* data_dir = std::getenv("SPDP_TEST_DATA_DIR");
    SPDP_CHECK(data_dir != nullptr);
    const auto data = spdp::read_spdp_data(
        std::string(data_dir) + "/" + filename
    );
    spdp::CpSatSolveOptions options;
    options.vehicle_count = vehicle_count;
    options.time_limit_seconds = 1200;
    options.workers = 0;
    const auto cp = spdp::solve_fixed_k_cp_sat(data, options);
    std::cerr << filename << " K=" << vehicle_count
              << " status=" << cp.status_name
              << " wall=" << cp.wall_time_seconds << '\n';
    SPDP_CHECK_EQ(cp.outcome, expected);
    if (require_adapter) {
        const auto graph = spdp::build_multigraph(data);
        const auto mapped = spdp::map_cp_incumbent_to_multigraph(
            data, graph, cp.routes
        );
        std::cerr << filename << " adapter=" << mapped.error_message << '\n';
        SPDP_CHECK(mapped.success);
    }
}

}  // namespace

SPDP_TEST(benchmark_d5_k7) {
    run_known_instance_case(
        "RecDep_day_D5.dat", 7,
        spdp::CpSolveOutcome::ProvenInfeasible, false
    );
}

SPDP_TEST(benchmark_d5_k8) {
    run_known_instance_case(
        "RecDep_day_D5.dat", 8,
        spdp::CpSolveOutcome::Feasible, true
    );
}

SPDP_TEST(benchmark_c19_k12) {
    run_known_instance_case(
        "RecDep_day_C19.dat", 12,
        spdp::CpSolveOutcome::ProvenInfeasible, false
    );
}

SPDP_TEST(benchmark_c19_k13) {
    const char* data_dir = std::getenv("SPDP_TEST_DATA_DIR");
    SPDP_CHECK(data_dir != nullptr);
    const auto data = spdp::read_spdp_data(
        std::string(data_dir) + "/RecDep_day_C19.dat"
    );
    spdp::CpSatSolveOptions options;
    options.vehicle_count = 13;
    options.time_limit_seconds = 1200;
    const auto cp = spdp::solve_fixed_k_cp_sat(data, options);
    SPDP_CHECK(cp.outcome == spdp::CpSolveOutcome::Feasible ||
               cp.outcome == spdp::CpSolveOutcome::Unknown);
    if (cp.outcome == spdp::CpSolveOutcome::Feasible) {
        const auto graph = spdp::build_multigraph(data);
        const auto mapped = spdp::map_cp_incumbent_to_multigraph(
            data, graph, cp.routes
        );
        SPDP_CHECK(mapped.success);
    }
}
