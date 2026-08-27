#include "CpSatSpdpSolver.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/sat_parameters.pb.h"

namespace spdp {
namespace {

using operations_research::Domain;
using operations_research::sat::BoolVar;
using operations_research::sat::CircuitConstraint;
using operations_research::sat::CpModelBuilder;
using operations_research::sat::CpSolverResponse;
using operations_research::sat::IntVar;
using operations_research::sat::LinearExpr;

constexpr double kIntegerTolerance = 1e-9;

struct IntegerData {
    std::int64_t pickup = 0;
    std::int64_t treatment = 0;
    std::int64_t delivery = 0;
    std::int64_t limit = 0;
    std::vector<std::vector<std::int64_t>> travel;
};

struct Action {
    CpActionKind kind = CpActionKind::Pickup;
    int request = -1;
    int location = -1;
    int type = -1;
    std::int64_t duration = 0;
};

struct ArcVariable {
    int tail = -1;
    int head = -1;
    BoolVar literal;
};

bool exact_integer(double value, std::int64_t& converted) {
    if (!std::isfinite(value)) {
        return false;
    }
    const double rounded = std::round(value);
    if (std::fabs(value - rounded) > kIntegerTolerance ||
        rounded < static_cast<double>(std::numeric_limits<std::int64_t>::min()) ||
        rounded > static_cast<double>(std::numeric_limits<std::int64_t>::max())) {
        return false;
    }
    converted = static_cast<std::int64_t>(rounded);
    return true;
}

bool convert_integer_data(
    const SPDPData& data,
    IntegerData& result,
    std::string& error
) {
    if (data.locations <= 0 ||
        static_cast<int>(data.time.size()) < data.locations) {
        error = "Time matrix has fewer rows than physical locations.";
        return false;
    }
    if (!exact_integer(data.time_pickup, result.pickup) ||
        !exact_integer(data.time_empty, result.treatment) ||
        !exact_integer(data.time_delivery, result.delivery) ||
        !exact_integer(data.time_limit, result.limit)) {
        error = "CP-SAT requires exact integer service times and route limit.";
        return false;
    }
    if (result.pickup < 0 || result.treatment < 0 || result.delivery < 0 ||
        result.limit < 0) {
        error = "CP-SAT time inputs must be nonnegative.";
        return false;
    }
    result.travel.assign(
        static_cast<std::size_t>(data.locations),
        std::vector<std::int64_t>(static_cast<std::size_t>(data.locations), 0)
    );
    for (int i = 0; i < data.locations; ++i) {
        if (static_cast<int>(data.time[static_cast<std::size_t>(i)].size()) <
            data.locations) {
            error = "Time matrix has fewer columns than physical locations.";
            return false;
        }
        for (int j = 0; j < data.locations; ++j) {
            std::int64_t value = 0;
            if (!exact_integer(
                    data.time[static_cast<std::size_t>(i)]
                             [static_cast<std::size_t>(j)],
                    value
                ) || value < 0) {
                error = "CP-SAT requires a nonnegative exact integer time matrix.";
                return false;
            }
            result.travel[static_cast<std::size_t>(i)]
                         [static_cast<std::size_t>(j)] = value;
        }
    }
    for (const Request& request : data.requests) {
        if (request.from_id < 0 || request.from_id >= data.locations ||
            request.to_id < 0 || request.to_id >= data.locations) {
            error = "Request references a location outside the time matrix.";
            return false;
        }
    }
    return true;
}

int action_index(CpActionKind kind, int request, int request_count) {
    switch (kind) {
    case CpActionKind::Pickup:
        return request;
    case CpActionKind::Treatment:
        return request_count + request;
    case CpActionKind::Delivery:
        return 2 * request_count + request;
    }
    return -1;
}

bool cor_40_forbidden(const Action& from, const Action& to) {
    return from.kind == to.kind && from.request > to.request &&
        from.location == to.location;
}

bool cor_41_forbidden(const Action& from, const Action& to) {
    return
        (from.kind == CpActionKind::Pickup &&
         to.kind == CpActionKind::Delivery &&
         from.location == to.location) ||
        (from.kind == CpActionKind::Treatment &&
         to.kind == CpActionKind::Pickup &&
         from.location == to.location) ||
        (from.kind == CpActionKind::Delivery &&
         to.kind == CpActionKind::Treatment &&
         from.location == to.location);
}

std::string action_name(const Action& action) {
    const char prefix = action.kind == CpActionKind::Pickup ? 'P' :
        (action.kind == CpActionKind::Treatment ? 'H' : 'D');
    return std::string(1, prefix) + std::to_string(action.request);
}

LinearExpr bool_sum(const std::vector<BoolVar>& vars) {
    LinearExpr expression;
    for (const BoolVar& var : vars) {
        expression += var;
    }
    return expression;
}

}  // namespace

CpFixedKSolveResult solve_fixed_k_cp_sat(
    const SPDPData& data,
    const CpSatSolveOptions& options
) {
    CpFixedKSolveResult result;
    if (options.vehicle_count <= 0 || options.workers < 0 ||
        options.random_seed < 0 || !std::isfinite(options.time_limit_seconds) ||
        options.time_limit_seconds < 0.0) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.status_name = "MODEL_INVALID";
        result.error_message = "Invalid CP-SAT solve options.";
        return result;
    }

    IntegerData integer_data;
    if (!convert_integer_data(data, integer_data, result.error_message)) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.status_name = "MODEL_INVALID";
        return result;
    }

    const int n = static_cast<int>(data.requests.size());
    const int k_count = options.vehicle_count;
    if (n == 0 || k_count > n) {
        result.outcome = CpSolveOutcome::ProvenInfeasible;
        result.status_name = "INFEASIBLE";
        return result;
    }
    const int action_count = 3 * n;

    std::vector<Action> actions(static_cast<std::size_t>(action_count));
    for (int i = 0; i < n; ++i) {
        const Request& request = data.requests[static_cast<std::size_t>(i)];
        actions[static_cast<std::size_t>(action_index(CpActionKind::Pickup, i, n))] =
            Action{CpActionKind::Pickup, i, request.from_id,
                   request.container_type, integer_data.pickup};
        actions[static_cast<std::size_t>(action_index(CpActionKind::Treatment, i, n))] =
            Action{CpActionKind::Treatment, i, request.to_id,
                   request.container_type, integer_data.treatment};
        actions[static_cast<std::size_t>(action_index(CpActionKind::Delivery, i, n))] =
            Action{CpActionKind::Delivery, i, request.from_id,
                   request.container_type, integer_data.delivery};
    }

    CpModelBuilder builder;
    std::vector<IntVar> vehicle;
    std::vector<IntVar> position;
    std::vector<IntVar> start_time;
    vehicle.reserve(static_cast<std::size_t>(action_count));
    position.reserve(static_cast<std::size_t>(action_count));
    start_time.reserve(static_cast<std::size_t>(action_count));
    for (const Action& action : actions) {
        vehicle.push_back(builder.NewIntVar(Domain(0, k_count - 1))
                              .WithName("vehicle_" + action_name(action)));
        position.push_back(builder.NewIntVar(Domain(1, action_count))
                               .WithName("position_" + action_name(action)));
        start_time.push_back(builder.NewIntVar(Domain(0, integer_data.limit))
                                 .WithName("time_" + action_name(action)));
    }

    std::vector<std::vector<BoolVar>> assigned(
        static_cast<std::size_t>(action_count),
        std::vector<BoolVar>(static_cast<std::size_t>(k_count))
    );
    for (int a = 0; a < action_count; ++a) {
        std::vector<BoolVar> row;
        for (int k = 0; k < k_count; ++k) {
            assigned[static_cast<std::size_t>(a)][static_cast<std::size_t>(k)] =
                builder.NewBoolVar().WithName(
                    "assigned_" + action_name(actions[static_cast<std::size_t>(a)]) +
                    "_" + std::to_string(k)
                );
            const BoolVar literal =
                assigned[static_cast<std::size_t>(a)][static_cast<std::size_t>(k)];
            builder.AddEquality(vehicle[static_cast<std::size_t>(a)], k)
                .OnlyEnforceIf(literal);
            row.push_back(literal);
        }
        builder.AddExactlyOne(row);
    }

    const int anchor_base = action_count;
    auto start_anchor = [anchor_base](int k) { return anchor_base + 2 * k; };
    auto end_anchor = [anchor_base](int k) { return anchor_base + 2 * k + 1; };
    CircuitConstraint circuit = builder.AddCircuitConstraint();
    std::vector<ArcVariable> arcs;
    std::vector<std::vector<BoolVar>> start_arcs(
        static_cast<std::size_t>(k_count));

    for (int k = 0; k < k_count; ++k) {
        for (int i = 0; i < n; ++i) {
            const int pickup = action_index(CpActionKind::Pickup, i, n);
            BoolVar literal = builder.NewBoolVar().WithName(
                "arc_S" + std::to_string(k) + "_P" + std::to_string(i)
            );
            circuit.AddArc(start_anchor(k), pickup, literal);
            arcs.push_back(ArcVariable{start_anchor(k), pickup, literal});
            start_arcs[static_cast<std::size_t>(k)].push_back(literal);
            builder.AddImplication(
                literal,
                assigned[static_cast<std::size_t>(pickup)][static_cast<std::size_t>(k)]
            );
            builder.AddEquality(position[static_cast<std::size_t>(pickup)], 1)
                .OnlyEnforceIf(literal);
            builder.AddEquality(start_time[static_cast<std::size_t>(pickup)], 0)
                .OnlyEnforceIf(literal);
        }
        for (int i = 0; i < n; ++i) {
            const int delivery = action_index(CpActionKind::Delivery, i, n);
            BoolVar literal = builder.NewBoolVar().WithName(
                "arc_D" + std::to_string(i) + "_Z" + std::to_string(k)
            );
            circuit.AddArc(delivery, end_anchor(k), literal);
            arcs.push_back(ArcVariable{delivery, end_anchor(k), literal});
            builder.AddImplication(
                literal,
                assigned[static_cast<std::size_t>(delivery)][static_cast<std::size_t>(k)]
            );
        }
        const int next_k = (k + 1) % k_count;
        circuit.AddArc(end_anchor(k), start_anchor(next_k), builder.TrueVar());
    }

    std::vector<ArcVariable> internal_arcs;
    for (int a = 0; a < action_count; ++a) {
        for (int b = 0; b < action_count; ++b) {
            if (a == b) {
                continue;
            }
            ++result.build_stats.candidate_action_arcs;
            const Action& from = actions[static_cast<std::size_t>(a)];
            const Action& to = actions[static_cast<std::size_t>(b)];
            if (from.kind == CpActionKind::Treatment &&
                to.kind == CpActionKind::Pickup &&
                from.request == to.request) {
                ++result.build_stats.core_precedence_pruned_arcs;
                continue;
            }
            if (cor_40_forbidden(from, to)) {
                ++result.build_stats.cor_40_pruned_arcs;
                continue;
            }
            if (cor_41_forbidden(from, to)) {
                ++result.build_stats.cor_41_pruned_arcs;
                continue;
            }
            BoolVar literal = builder.NewBoolVar().WithName(
                "arc_" + action_name(from) + "_" + action_name(to)
            );
            circuit.AddArc(a, b, literal);
            ArcVariable arc{a, b, literal};
            arcs.push_back(arc);
            internal_arcs.push_back(arc);
            ++result.build_stats.created_action_arcs;
            builder.AddEquality(vehicle[static_cast<std::size_t>(b)],
                                vehicle[static_cast<std::size_t>(a)])
                .OnlyEnforceIf(literal);
            builder.AddEquality(position[static_cast<std::size_t>(b)],
                                position[static_cast<std::size_t>(a)] + 1)
                .OnlyEnforceIf(literal);
            const std::int64_t transition =
                from.duration + integer_data.travel[static_cast<std::size_t>(from.location)]
                                                   [static_cast<std::size_t>(to.location)];
            builder.AddEquality(start_time[static_cast<std::size_t>(b)],
                                start_time[static_cast<std::size_t>(a)] + transition)
                .OnlyEnforceIf(literal);
        }
    }

    std::vector<IntVar> completion;
    completion.reserve(static_cast<std::size_t>(k_count));
    for (int k = 0; k < k_count; ++k) {
        completion.push_back(builder.NewIntVar(Domain(0, integer_data.limit))
                                 .WithName("completion_" + std::to_string(k)));
        for (const ArcVariable& arc : arcs) {
            if (arc.head != end_anchor(k)) {
                continue;
            }
            builder.AddEquality(
                completion[static_cast<std::size_t>(k)],
                start_time[static_cast<std::size_t>(arc.tail)] +
                    actions[static_cast<std::size_t>(arc.tail)].duration
            ).OnlyEnforceIf(arc.literal);
        }
    }

    for (int i = 0; i < n; ++i) {
        const int pickup = action_index(CpActionKind::Pickup, i, n);
        const int treatment = action_index(CpActionKind::Treatment, i, n);
        builder.AddEquality(vehicle[static_cast<std::size_t>(pickup)],
                            vehicle[static_cast<std::size_t>(treatment)]);
        builder.AddLessOrEqual(position[static_cast<std::size_t>(pickup)] + 1,
                               position[static_cast<std::size_t>(treatment)]);
    }

    std::set<int> types;
    for (const Request& request : data.requests) {
        types.insert(request.container_type);
    }
    for (int k = 0; k < k_count; ++k) {
        auto load = builder.AddReservoirConstraint(0, 2);
        for (int i = 0; i < n; ++i) {
            const int pickup = action_index(CpActionKind::Pickup, i, n);
            const int delivery = action_index(CpActionKind::Delivery, i, n);
            load.AddOptionalEvent(position[static_cast<std::size_t>(pickup)], 1,
                assigned[static_cast<std::size_t>(pickup)][static_cast<std::size_t>(k)]);
            load.AddOptionalEvent(position[static_cast<std::size_t>(delivery)], -1,
                assigned[static_cast<std::size_t>(delivery)][static_cast<std::size_t>(k)]);
        }
        for (int type : types) {
            auto empty_inventory = builder.AddReservoirConstraint(0, 2);
            for (int i = 0; i < n; ++i) {
                if (data.requests[static_cast<std::size_t>(i)].container_type != type) {
                    continue;
                }
                const int treatment = action_index(CpActionKind::Treatment, i, n);
                const int delivery = action_index(CpActionKind::Delivery, i, n);
                empty_inventory.AddOptionalEvent(
                    position[static_cast<std::size_t>(treatment)], 1,
                    assigned[static_cast<std::size_t>(treatment)][static_cast<std::size_t>(k)]
                );
                empty_inventory.AddOptionalEvent(
                    position[static_cast<std::size_t>(delivery)], -1,
                    assigned[static_cast<std::size_t>(delivery)][static_cast<std::size_t>(k)]
                );
            }
        }
        if (options.redundant.full_skip_reservoir) {
            auto full_inventory = builder.AddReservoirConstraint(0, 2);
            for (int i = 0; i < n; ++i) {
                const int pickup = action_index(CpActionKind::Pickup, i, n);
                const int treatment = action_index(CpActionKind::Treatment, i, n);
                full_inventory.AddOptionalEvent(
                    position[static_cast<std::size_t>(pickup)], 1,
                    assigned[static_cast<std::size_t>(pickup)][static_cast<std::size_t>(k)]
                );
                full_inventory.AddOptionalEvent(
                    position[static_cast<std::size_t>(treatment)], -1,
                    assigned[static_cast<std::size_t>(treatment)][static_cast<std::size_t>(k)]
                );
            }
        }
    }

    if (options.redundant.terminal_balance) {
        for (int k = 0; k < k_count; ++k) {
            LinearExpr pickups;
            LinearExpr deliveries;
            for (int i = 0; i < n; ++i) {
                pickups += assigned[static_cast<std::size_t>(
                    action_index(CpActionKind::Pickup, i, n))][static_cast<std::size_t>(k)];
                deliveries += assigned[static_cast<std::size_t>(
                    action_index(CpActionKind::Delivery, i, n))][static_cast<std::size_t>(k)];
            }
            builder.AddEquality(pickups, deliveries);
            for (int type : types) {
                LinearExpr treatments;
                LinearExpr typed_deliveries;
                for (int i = 0; i < n; ++i) {
                    if (data.requests[static_cast<std::size_t>(i)].container_type != type) {
                        continue;
                    }
                    treatments += assigned[static_cast<std::size_t>(
                        action_index(CpActionKind::Treatment, i, n))][static_cast<std::size_t>(k)];
                    typed_deliveries += assigned[static_cast<std::size_t>(
                        action_index(CpActionKind::Delivery, i, n))][static_cast<std::size_t>(k)];
                }
                builder.AddEquality(treatments, typed_deliveries);
            }
        }
    }

    if (options.redundant.container_workload) {
        for (int k = 0; k < k_count; ++k) {
            LinearExpr service_twice;
            LinearExpr workload;
            for (int i = 0; i < n; ++i) {
                const Request& request = data.requests[static_cast<std::size_t>(i)];
                const int pickup = action_index(CpActionKind::Pickup, i, n);
                const int treatment = action_index(CpActionKind::Treatment, i, n);
                const int delivery = action_index(CpActionKind::Delivery, i, n);
                service_twice += 2 * integer_data.pickup *
                    assigned[static_cast<std::size_t>(pickup)][static_cast<std::size_t>(k)];
                service_twice += 2 * integer_data.treatment *
                    assigned[static_cast<std::size_t>(treatment)][static_cast<std::size_t>(k)];
                service_twice += 2 * integer_data.delivery *
                    assigned[static_cast<std::size_t>(delivery)][static_cast<std::size_t>(k)];
                workload += integer_data.travel[static_cast<std::size_t>(request.from_id)]
                                                     [static_cast<std::size_t>(request.to_id)] *
                    assigned[static_cast<std::size_t>(pickup)][static_cast<std::size_t>(k)];
            }
            std::vector<std::vector<BoolVar>> matching(
                static_cast<std::size_t>(n),
                std::vector<BoolVar>(static_cast<std::size_t>(n))
            );
            for (int i = 0; i < n; ++i) {
                std::vector<BoolVar> outgoing;
                for (int j = 0; j < n; ++j) {
                    if (data.requests[static_cast<std::size_t>(i)].container_type !=
                        data.requests[static_cast<std::size_t>(j)].container_type) {
                        continue;
                    }
                    matching[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                        builder.NewBoolVar().WithName(
                            "lambda_" + std::to_string(i) + "_" +
                            std::to_string(j) + "_" + std::to_string(k)
                        );
                    const BoolVar lambda = matching[static_cast<std::size_t>(i)]
                                                  [static_cast<std::size_t>(j)];
                    outgoing.push_back(lambda);
                    workload += integer_data.travel[static_cast<std::size_t>(
                        data.requests[static_cast<std::size_t>(i)].to_id)]
                        [static_cast<std::size_t>(
                            data.requests[static_cast<std::size_t>(j)].from_id)] * lambda;
                }
                builder.AddEquality(
                    bool_sum(outgoing),
                    assigned[static_cast<std::size_t>(
                        action_index(CpActionKind::Pickup, i, n))][static_cast<std::size_t>(k)]
                );
            }
            for (int j = 0; j < n; ++j) {
                std::vector<BoolVar> incoming;
                for (int i = 0; i < n; ++i) {
                    if (data.requests[static_cast<std::size_t>(i)].container_type ==
                        data.requests[static_cast<std::size_t>(j)].container_type) {
                        incoming.push_back(matching[static_cast<std::size_t>(i)]
                                                  [static_cast<std::size_t>(j)]);
                    }
                }
                builder.AddEquality(
                    bool_sum(incoming),
                    assigned[static_cast<std::size_t>(
                        action_index(CpActionKind::Delivery, j, n))][static_cast<std::size_t>(k)]
                );
            }
            builder.AddLessOrEqual(service_twice + workload, 2 * integer_data.limit);
        }
    }

    if (options.redundant.aggregate_duration) {
        LinearExpr total_duration =
            n * (integer_data.pickup + integer_data.treatment + integer_data.delivery);
        for (const ArcVariable& arc : internal_arcs) {
            const Action& from = actions[static_cast<std::size_t>(arc.tail)];
            const Action& to = actions[static_cast<std::size_t>(arc.head)];
            total_duration +=
                integer_data.travel[static_cast<std::size_t>(from.location)]
                                   [static_cast<std::size_t>(to.location)] * arc.literal;
        }
        builder.AddLessOrEqual(total_duration, k_count * integer_data.limit);
    }

    if (options.symmetry.first_pickup_vehicle_ordering) {
        for (int k = 0; k + 1 < k_count; ++k) {
            LinearExpr current;
            LinearExpr next;
            for (int i = 0; i < n; ++i) {
                current += (i + 1) * start_arcs[static_cast<std::size_t>(k)]
                                             [static_cast<std::size_t>(i)];
                next += (i + 1) * start_arcs[static_cast<std::size_t>(k + 1)]
                                          [static_cast<std::size_t>(i)];
            }
            builder.AddLessOrEqual(current + 1, next);
        }
    }

    if (options.symmetry.cor_43_identical_pickup_time_ordering) {
        std::map<std::tuple<int, int, int>, std::vector<int>> groups;
        for (int i = 0; i < n; ++i) {
            const Request& request = data.requests[static_cast<std::size_t>(i)];
            groups[{request.from_id, request.to_id, request.container_type}]
                .push_back(i);
        }
        for (const auto& entry : groups) {
            const std::vector<int>& requests = entry.second;
            for (std::size_t left = 0; left < requests.size(); ++left) {
                for (std::size_t right = left + 1; right < requests.size(); ++right) {
                    builder.AddLessOrEqual(
                        start_time[static_cast<std::size_t>(action_index(
                            CpActionKind::Pickup, requests[left], n))],
                        start_time[static_cast<std::size_t>(action_index(
                            CpActionKind::Pickup, requests[right], n))]
                    );
                    ++result.build_stats.cor_43_constraint_count;
                }
            }
        }
    }

    operations_research::sat::SatParameters parameters;
    if (options.time_limit_seconds > 0.0) {
        parameters.set_max_time_in_seconds(options.time_limit_seconds);
    }
    if (options.workers > 0) {
        parameters.set_num_search_workers(options.workers);
    }
    parameters.set_random_seed(options.random_seed);
    parameters.set_log_search_progress(options.log_search_progress);
    operations_research::sat::Model model;
    model.Add(operations_research::sat::NewSatParameters(parameters));
    const CpSolverResponse response =
        operations_research::sat::SolveCpModel(builder.Build(), &model);

    result.raw_status = static_cast<int>(response.status());
    result.status_name = operations_research::sat::CpSolverStatus_Name(response.status());
    result.wall_time_seconds = response.wall_time();
    result.conflicts = response.num_conflicts();
    result.branches = response.num_branches();
    if (response.status() == operations_research::sat::CpSolverStatus::INFEASIBLE) {
        result.outcome = CpSolveOutcome::ProvenInfeasible;
        return result;
    }
    if (response.status() == operations_research::sat::CpSolverStatus::MODEL_INVALID) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.error_message = response.solution_info();
        return result;
    }
    if (response.status() != operations_research::sat::CpSolverStatus::FEASIBLE &&
        response.status() != operations_research::sat::CpSolverStatus::OPTIMAL) {
        result.outcome = CpSolveOutcome::Unknown;
        return result;
    }

    result.outcome = CpSolveOutcome::Feasible;
    std::map<int, int> successor;
    for (const ArcVariable& arc : arcs) {
        if (operations_research::sat::SolutionBooleanValue(response, arc.literal)) {
            successor[arc.tail] = arc.head;
        }
    }
    std::vector<bool> visited(static_cast<std::size_t>(action_count), false);
    for (int k = 0; k < k_count; ++k) {
        CpActionRoute route;
        route.vehicle_index = k;
        int current = start_anchor(k);
        for (int step = 0; step <= action_count; ++step) {
            const auto found = successor.find(current);
            if (found == successor.end()) {
                result.outcome = CpSolveOutcome::ModelInvalid;
                result.error_message = "Selected CP circuit cannot be decoded.";
                result.routes.clear();
                return result;
            }
            current = found->second;
            if (current == end_anchor(k)) {
                break;
            }
            if (current < 0 || current >= action_count ||
                visited[static_cast<std::size_t>(current)]) {
                result.outcome = CpSolveOutcome::ModelInvalid;
                result.error_message = "Selected CP circuit repeats or leaves an action route.";
                result.routes.clear();
                return result;
            }
            visited[static_cast<std::size_t>(current)] = true;
            const Action& action = actions[static_cast<std::size_t>(current)];
            route.actions.push_back(CpActionVisit{
                action.kind,
                action.request,
                static_cast<int>(operations_research::sat::SolutionIntegerValue(
                    response, position[static_cast<std::size_t>(current)])),
                operations_research::sat::SolutionIntegerValue(
                    response, start_time[static_cast<std::size_t>(current)])
            });
        }
        route.completion_time = operations_research::sat::SolutionIntegerValue(
            response, completion[static_cast<std::size_t>(k)]
        );
        result.routes.push_back(std::move(route));
    }
    if (std::count(visited.begin(), visited.end(), true) != action_count) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.error_message = "Selected CP circuit does not cover every action exactly once.";
        result.routes.clear();
    }
    return result;
}

}  // namespace spdp
