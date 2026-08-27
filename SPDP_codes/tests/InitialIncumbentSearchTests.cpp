#include "TestSupport.h"

#include "InitialIncumbent.h"

#include <filesystem>

namespace {

spdp::InitialIncumbentSolveResult result(
    int k,
    spdp::FixedKSolveOutcome outcome
) {
    spdp::InitialIncumbentSolveResult value;
    value.vehicle_count = k;
    value.outcome = outcome;
    value.has_feasible_solution = outcome == spdp::FixedKSolveOutcome::Feasible;
    return value;
}

}  // namespace

SPDP_TEST(initial_incumbent_search) {
    spdp::InitialIncumbentSearchOptions options;
    SPDP_CHECK_EQ(
        options.backend,
        spdp::InitialIncumbentBackend::TwoIndexMilp
    );
    options.initial_vehicle_count = 8;
    options.max_k_increments = 1;
    options.timeout_action = spdp::InitialIncumbentTimeoutAction::Advance;

    int call = 0;
    const auto unknown_then_feasible =
        spdp::run_iterative_initial_incumbent_search(
            options,
            20,
            [&call](int k) {
                ++call;
                return call == 1
                    ? result(k, spdp::FixedKSolveOutcome::Unknown)
                    : result(k, spdp::FixedKSolveOutcome::Feasible);
            }
        );
    SPDP_CHECK_EQ(unknown_then_feasible.certified_k, 8);
    SPDP_CHECK_EQ(unknown_then_feasible.incumbent_k, 9);

    options.initial_vehicle_count = 7;
    options.max_k_increments = 0;
    options.timeout_action = spdp::InitialIncumbentTimeoutAction::Stop;
    const auto infeasible = spdp::run_iterative_initial_incumbent_search(
        options,
        20,
        [](int k) {
            return result(k, spdp::FixedKSolveOutcome::ProvenInfeasible);
        }
    );
    SPDP_CHECK_EQ(infeasible.certified_k, 8);

    options.initial_vehicle_count = 3;
    options.max_k_increments = 2;
    int invalid_calls = 0;
    const auto invalid = spdp::run_iterative_initial_incumbent_search(
        options,
        10,
        [&invalid_calls](int k) {
            ++invalid_calls;
            return result(k, spdp::FixedKSolveOutcome::ModelInvalid);
        }
    );
    SPDP_CHECK_EQ(invalid_calls, 1);
    SPDP_CHECK_EQ(invalid.certified_k, 3);
}

SPDP_TEST(initial_incumbent_backend) {
    spdp::SPDPData data;
    data.fixed_vehicle_cost = 500;
    data.time_pickup = 2;
    data.time_empty = 3;
    data.time_delivery = 2;
    data.time_limit = 100;
    data.locations = 3;
    data.requests = {
        spdp::Request{0, 2, 4},
        spdp::Request{1, 2, 4},
    };
    data.time = {{0, 2, 5, 0}, {2, 0, 4, 0}, {5, 4, 0, 0}, {0, 0, 0, 0}};
    data.distance = data.time;
    const auto graph = spdp::build_multigraph(data);

    spdp::InitialIncumbentSearchOptions options;
    options.backend = spdp::InitialIncumbentBackend::CpSat;
    options.initial_vehicle_count = 1;
    options.max_k_increments = 0;
    options.per_attempt_time_limit = 5;
    const std::filesystem::path log_base =
        std::filesystem::temp_directory_path() /
        "RecDep_day_test_initial_incumbent_cp_sat.log";
    const std::filesystem::path expected_log =
        std::filesystem::temp_directory_path() /
        "RecDep_day_test_initial_incumbent_k1_cp_sat.log";
    std::filesystem::remove(log_base);
    std::filesystem::remove(expected_log);
    options.cp_sat_log_base_path = log_base.string();
    const auto search = spdp::solve_iterative_duration_initial_incumbent(
        data, graph, options
    );
    SPDP_CHECK(search.has_feasible_solution);
    SPDP_CHECK_EQ(search.incumbent_k, 1);
    SPDP_CHECK_EQ(
        search.incumbent_result.outcome,
        spdp::FixedKSolveOutcome::Feasible
    );
    SPDP_CHECK_EQ(search.incumbent_result.adapter_status, std::string("passed"));
    SPDP_CHECK(std::filesystem::exists(expected_log));
    std::filesystem::remove(expected_log);
}
