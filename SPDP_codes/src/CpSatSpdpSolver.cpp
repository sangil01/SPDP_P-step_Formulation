#include "CpSatSpdpSolver.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "ortools/util/logging.h"

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
constexpr double kObjectiveTolerance = 1e-9;

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

bool compute_threshold_horizon(
    std::int64_t limit,
    double factor,
    std::int64_t& horizon,
    std::string& error
) {
    if (limit <= 0 || !std::isfinite(factor) || factor <= 1.0) {
        error =
            "Threshold optimization requires a positive route limit and "
            "a finite horizon factor greater than one.";
        return false;
    }
    const long double scaled =
        static_cast<long double>(limit) * static_cast<long double>(factor);
    if (!std::isfinite(scaled) ||
        scaled > static_cast<long double>(std::numeric_limits<std::int64_t>::max())) {
        error = "Threshold horizon exceeds the CP-SAT integer range.";
        return false;
    }
    horizon = static_cast<std::int64_t>(std::ceil(scaled));
    if (horizon <= limit) {
        error = "Threshold horizon must be strictly greater than the route limit.";
        return false;
    }
    return true;
}

bool objective_bound_excludes_threshold(double bound, std::int64_t limit) {
    return std::isfinite(bound) &&
        bound > static_cast<double>(limit) + kObjectiveTolerance;
}

}  // namespace

const char* cp_solve_mode_name(CpSolveMode mode) {
    switch (mode) {
    case CpSolveMode::Satisfaction:
        return "satisfaction";
    case CpSolveMode::ThresholdOptimization:
        return "threshold-optimization";
    }
    return "unknown";
}

const char* cp_graph_mode_name(CpGraphMode mode) {
    switch (mode) {
    case CpGraphMode::OriginalGraph:
        return "original-graph";
    case CpGraphMode::Multigraph:
        return "multigraph";
    }
    return "unknown";
}

CpSolveOutcome classify_cp_unknown_outcome(
    double configured_time_limit_seconds,
    double solver_wall_time_seconds
) {
    if (!(configured_time_limit_seconds > 0.0) ||
        !std::isfinite(configured_time_limit_seconds) ||
        !std::isfinite(solver_wall_time_seconds)) {
        return CpSolveOutcome::EarlyUnknown;
    }

    // CP-SAT can report a wall time a few milliseconds below its configured
    // limit.  Cap the tolerance at 0.1 s so substantially early returns such
    // as 87.47/1200 s can never be classified as timeouts.
    const double tolerance = std::min(
        0.1,
        std::max(0.001, configured_time_limit_seconds * 1e-6)
    );
    return solver_wall_time_seconds + tolerance >= configured_time_limit_seconds
        ? CpSolveOutcome::TimedOutUnknown
        : CpSolveOutcome::EarlyUnknown;
}

static CpFixedKSolveResult solve_fixed_k_original_graph_cp_sat(
    const SPDPData& data,
    const CpSatSolveOptions& options
) {
    CpFixedKSolveResult result;
    result.graph_mode = CpGraphMode::OriginalGraph;
    result.solve_mode = options.solve_mode;
    result.configured_time_limit_seconds = options.time_limit_seconds;
    result.threshold_horizon_factor =
        options.solve_mode == CpSolveMode::ThresholdOptimization
            ? options.threshold_horizon_factor
            : 0.0;
    const bool write_log_file = !options.log_file_path.empty();
    std::ofstream log_file;
    if (write_log_file) {
        log_file.open(options.log_file_path, std::ios::trunc);
        if (!log_file) {
            result.outcome = CpSolveOutcome::ModelInvalid;
            result.status_name = "MODEL_INVALID";
            result.termination_name = "model-invalid";
            result.error_message =
                "Failed to create CP-SAT log file: " + options.log_file_path;
            return result;
        }
        log_file << "Starting CP-SAT fixed-K incumbent attempt\n";
        log_file.flush();
    }
    auto append_early_log = [&log_file, write_log_file](
                                const std::string& status,
                                const std::string& message) {
        if (!write_log_file) {
            return;
        }
        log_file << status << ": " << message << '\n';
        log_file.flush();
    };
    if (options.vehicle_count <= 0 || options.workers < 0 ||
        !std::isfinite(options.time_limit_seconds) ||
        options.time_limit_seconds < 0.0) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.status_name = "MODEL_INVALID";
        result.termination_name = "model-invalid";
        result.error_message = "Invalid CP-SAT solve options.";
        append_early_log(result.status_name, result.error_message);
        return result;
    }

    IntegerData integer_data;
    if (!convert_integer_data(data, integer_data, result.error_message)) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.status_name = "MODEL_INVALID";
        result.termination_name = "model-invalid";
        append_early_log(result.status_name, result.error_message);
        return result;
    }

    std::int64_t route_horizon = integer_data.limit;
    if (options.solve_mode == CpSolveMode::ThresholdOptimization &&
        !compute_threshold_horizon(
            integer_data.limit,
            options.threshold_horizon_factor,
            route_horizon,
            result.error_message
        )) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.status_name = "MODEL_INVALID";
        result.termination_name = "model-invalid";
        append_early_log(result.status_name, result.error_message);
        return result;
    }
    result.model_horizon = route_horizon;

    const std::int64_t largest_rhs_multiplier =
        std::max<std::int64_t>(2, options.vehicle_count);
    if (route_horizon >
        std::numeric_limits<std::int64_t>::max() / largest_rhs_multiplier) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.status_name = "MODEL_INVALID";
        result.termination_name = "model-invalid";
        result.error_message = "Active CP horizon makes a redundant bound overflow.";
        append_early_log(result.status_name, result.error_message);
        return result;
    }

    const int n = static_cast<int>(data.requests.size());
    const int k_count = options.vehicle_count;
    if (n == 0 || k_count > n) {
        result.outcome = CpSolveOutcome::ProvenInfeasible;
        result.status_name = "INFEASIBLE";
        result.termination_name = "proven-infeasible";
        append_early_log(
            result.status_name,
            "Exact vehicle count exceeds the number of requests."
        );
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
    for (const Action& action : actions) {
        if (action.duration > route_horizon) {
            result.outcome = CpSolveOutcome::ProvenInfeasible;
            result.status_name = "INFEASIBLE";
            result.termination_name =
                options.solve_mode == CpSolveMode::ThresholdOptimization
                    ? "restricted-horizon-infeasible"
                    : "proven-infeasible";
            append_early_log(
                result.status_name,
                "A mandatory action exceeds the active route horizon."
            );
            return result;
        }
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
        start_time.push_back(builder.NewIntVar(Domain(
                                 0, route_horizon - action.duration))
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
        completion.push_back(builder.NewIntVar(Domain(0, route_horizon))
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

    std::optional<IntVar> makespan;
    if (options.solve_mode == CpSolveMode::ThresholdOptimization) {
        makespan.emplace(
            builder.NewIntVar(Domain(0, route_horizon)).WithName("makespan")
        );
        builder.AddMaxEquality(*makespan, completion);
        builder.Minimize(*makespan);
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
            builder.AddLessOrEqual(service_twice + workload, 2 * route_horizon);
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
        builder.AddLessOrEqual(total_duration, k_count * route_horizon);
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
    parameters.set_log_search_progress(write_log_file);
    parameters.set_log_to_stdout(false);
    operations_research::sat::Model model;
    model.Add(operations_research::sat::NewSatParameters(parameters));
    if (write_log_file) {
        operations_research::SolverLogger* logger =
            model.GetOrCreate<operations_research::SolverLogger>();
        logger->EnableLogging(true);
        logger->SetLogToStdOut(false);
        logger->AddInfoLoggingCallback([&log_file](const std::string& message) {
            log_file << message;
            if (message.empty() || message.back() != '\n') {
                log_file << '\n';
            }
            log_file.flush();
        });
    }

    std::mutex callback_mutex;
    std::optional<CpSolverResponse> accepted_threshold_response;
    bool stopped_by_bound_callback = false;
    double callback_proving_bound = 0.0;
    if (options.solve_mode == CpSolveMode::ThresholdOptimization) {
        model.Add(operations_research::sat::NewFeasibleSolutionObserver(
            [&](const CpSolverResponse& candidate) {
                if (candidate.objective_value() >
                    static_cast<double>(integer_data.limit) +
                        kObjectiveTolerance) {
                    return;
                }
                bool should_stop = false;
                {
                    std::lock_guard<std::mutex> lock(callback_mutex);
                    if (!accepted_threshold_response.has_value()) {
                        accepted_threshold_response = candidate;
                        should_stop = true;
                    }
                }
                if (should_stop) {
                    operations_research::sat::StopSearch(&model);
                }
            }
        ));
        model.Add(operations_research::sat::NewBestBoundCallback(
            [&](double bound) {
                if (!objective_bound_excludes_threshold(
                        bound, integer_data.limit)) {
                    return;
                }
                bool should_stop = false;
                {
                    std::lock_guard<std::mutex> lock(callback_mutex);
                    if (!stopped_by_bound_callback ||
                        bound > callback_proving_bound) {
                        callback_proving_bound = bound;
                    }
                    if (!stopped_by_bound_callback) {
                        stopped_by_bound_callback = true;
                        should_stop = true;
                    }
                }
                if (should_stop) {
                    operations_research::sat::StopSearch(&model);
                }
            }
        ));
    }

    const CpSolverResponse response =
        operations_research::sat::SolveCpModel(builder.Build(), &model);

    std::optional<CpSolverResponse> saved_threshold_response;
    {
        std::lock_guard<std::mutex> lock(callback_mutex);
        saved_threshold_response = accepted_threshold_response;
        result.stopped_by_feasible_observer =
            accepted_threshold_response.has_value();
        result.stopped_by_bound_callback = stopped_by_bound_callback;
        if (stopped_by_bound_callback) {
            result.has_objective_bound = true;
            result.best_objective_bound = callback_proving_bound;
        }
    }

    result.raw_status = static_cast<int>(response.status());
    result.status_name = operations_research::sat::CpSolverStatus_Name(response.status());
    result.solution_info = response.solution_info();
    result.wall_time_seconds = response.wall_time();
    result.conflicts = response.num_conflicts();
    result.branches = response.num_branches();
    const bool raw_has_solution =
        response.status() == operations_research::sat::CpSolverStatus::FEASIBLE ||
        response.status() == operations_research::sat::CpSolverStatus::OPTIMAL;
    if (options.solve_mode == CpSolveMode::ThresholdOptimization) {
        if (raw_has_solution) {
            result.has_objective_value = true;
            result.objective_value = response.objective_value();
        }
        if (std::isfinite(response.best_objective_bound())) {
            if (!result.has_objective_bound ||
                response.best_objective_bound() > result.best_objective_bound) {
                result.best_objective_bound = response.best_objective_bound();
            }
            result.has_objective_bound = true;
        }
        if (saved_threshold_response.has_value()) {
            result.has_objective_value = true;
            result.objective_value =
                saved_threshold_response->objective_value();
        }
    }

    auto append_final_log = [&]() {
        if (!write_log_file) {
            return;
        }
        log_file << "[spdp-cp-wrapper] graph_mode=original-graph mode="
                 << cp_solve_mode_name(result.solve_mode)
                 << " horizon_factor=" << result.threshold_horizon_factor
                 << " model_horizon=" << result.model_horizon
                 << " termination=" << result.termination_name
                 << " has_objective_value="
                 << (result.has_objective_value ? 1 : 0)
                 << " objective_value=" << result.objective_value
                 << " has_objective_bound="
                 << (result.has_objective_bound ? 1 : 0)
                 << " best_objective_bound=" << result.best_objective_bound
                 << " stopped_by_feasible_observer="
                 << (result.stopped_by_feasible_observer ? 1 : 0)
                 << " stopped_by_bound_callback="
                 << (result.stopped_by_bound_callback ? 1 : 0)
                 << '\n';
        log_file.flush();
    };

    if (response.status() == operations_research::sat::CpSolverStatus::MODEL_INVALID) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.termination_name = "model-invalid";
        result.error_message = result.solution_info;
        append_final_log();
        return result;
    }

    const CpSolverResponse* solution_response = nullptr;
    if (options.solve_mode == CpSolveMode::Satisfaction) {
        if (response.status() == operations_research::sat::CpSolverStatus::INFEASIBLE) {
            result.outcome = CpSolveOutcome::ProvenInfeasible;
            result.termination_name = "proven-infeasible";
            append_final_log();
            return result;
        }
        if (raw_has_solution) {
            result.outcome = CpSolveOutcome::Feasible;
            result.termination_name = "feasible";
            solution_response = &response;
        }
    } else {
        if (saved_threshold_response.has_value()) {
            result.outcome = CpSolveOutcome::Feasible;
            result.termination_name = "threshold-feasible-witness";
            solution_response = &*saved_threshold_response;
        } else if (raw_has_solution &&
                   response.objective_value() <=
                       static_cast<double>(integer_data.limit) +
                           kObjectiveTolerance) {
            result.outcome = CpSolveOutcome::Feasible;
            result.termination_name = "threshold-feasible-witness";
            solution_response = &response;
        } else if (result.stopped_by_bound_callback ||
                   (result.has_objective_bound &&
                    objective_bound_excludes_threshold(
                        result.best_objective_bound, integer_data.limit)) ||
                   (response.status() ==
                        operations_research::sat::CpSolverStatus::OPTIMAL &&
                    raw_has_solution &&
                    response.objective_value() >
                        static_cast<double>(integer_data.limit) +
                            kObjectiveTolerance)) {
            result.outcome = CpSolveOutcome::ProvenInfeasible;
            result.termination_name = "objective-bound-infeasible";
            append_final_log();
            return result;
        } else if (response.status() ==
                   operations_research::sat::CpSolverStatus::INFEASIBLE) {
            result.outcome = CpSolveOutcome::ProvenInfeasible;
            result.termination_name = "restricted-horizon-infeasible";
            append_final_log();
            return result;
        }
    }

    if (solution_response == nullptr) {
        result.outcome = classify_cp_unknown_outcome(
            result.configured_time_limit_seconds,
            result.wall_time_seconds
        );
        result.hit_time_limit = result.outcome == CpSolveOutcome::TimedOutUnknown;
        result.early_unknown = result.outcome == CpSolveOutcome::EarlyUnknown;
        result.termination_name = result.hit_time_limit
            ? "timed-out-unknown" : "early-unknown";
        append_final_log();
        return result;
    }

    std::map<int, int> successor;
    for (const ArcVariable& arc : arcs) {
        if (operations_research::sat::SolutionBooleanValue(
                *solution_response, arc.literal)) {
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
                result.termination_name = "model-invalid";
                result.error_message = "Selected CP circuit cannot be decoded.";
                result.routes.clear();
                append_final_log();
                return result;
            }
            current = found->second;
            if (current == end_anchor(k)) {
                break;
            }
            if (current < 0 || current >= action_count ||
                visited[static_cast<std::size_t>(current)]) {
                result.outcome = CpSolveOutcome::ModelInvalid;
                result.termination_name = "model-invalid";
                result.error_message = "Selected CP circuit repeats or leaves an action route.";
                result.routes.clear();
                append_final_log();
                return result;
            }
            visited[static_cast<std::size_t>(current)] = true;
            const Action& action = actions[static_cast<std::size_t>(current)];
            route.actions.push_back(CpActionVisit{
                action.kind,
                action.request,
                static_cast<int>(operations_research::sat::SolutionIntegerValue(
                    *solution_response,
                    position[static_cast<std::size_t>(current)])),
                operations_research::sat::SolutionIntegerValue(
                    *solution_response,
                    start_time[static_cast<std::size_t>(current)])
            });
        }
        route.completion_time = operations_research::sat::SolutionIntegerValue(
            *solution_response, completion[static_cast<std::size_t>(k)]
        );
        result.routes.push_back(std::move(route));
    }
    if (std::count(visited.begin(), visited.end(), true) != action_count) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.termination_name = "model-invalid";
        result.error_message = "Selected CP circuit does not cover every action exactly once.";
        result.routes.clear();
    }
    append_final_log();
    return result;
}

static CpFixedKSolveResult solve_fixed_k_multigraph_cp_sat(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const CpSatSolveOptions& options
) {
    struct MultigraphArcVariable {
        int tail = -1;
        int head = -1;
        int original_edge_id = -1;
        int vehicle_index = -1;
        BoolVar literal;
    };
    struct StateFlowTerms {
        std::vector<BoolVar> incoming;
        std::vector<BoolVar> outgoing;
    };

    CpFixedKSolveResult result;
    result.graph_mode = CpGraphMode::Multigraph;
    result.solve_mode = options.solve_mode;
    result.configured_time_limit_seconds = options.time_limit_seconds;
    result.threshold_horizon_factor =
        options.solve_mode == CpSolveMode::ThresholdOptimization
            ? options.threshold_horizon_factor
            : 0.0;
    result.build_stats.input_multigraph_edges =
        static_cast<std::int64_t>(graph.number_of_edges());
    result.build_stats.full_skip_state_embedded = true;

    const bool write_log_file = !options.log_file_path.empty();
    std::ofstream log_file;
    if (write_log_file) {
        log_file.open(options.log_file_path, std::ios::trunc);
        if (!log_file) {
            result.outcome = CpSolveOutcome::ModelInvalid;
            result.status_name = "MODEL_INVALID";
            result.termination_name = "model-invalid";
            result.error_message =
                "Failed to create CP-SAT log file: " + options.log_file_path;
            return result;
        }
        log_file << "Starting CP-SAT fixed-K multigraph incumbent attempt\n";
        log_file.flush();
    }
    auto fail_model = [&](const std::string& message) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.status_name = "MODEL_INVALID";
        result.termination_name = "model-invalid";
        result.error_message = message;
        if (write_log_file) {
            log_file << result.status_name << ": " << message << '\n';
            log_file.flush();
        }
        return result;
    };

    if (options.vehicle_count <= 0 || options.workers < 0 ||
        !std::isfinite(options.time_limit_seconds) ||
        options.time_limit_seconds < 0.0) {
        return fail_model("Invalid CP-SAT solve options.");
    }

    IntegerData integer_data;
    std::string integer_error;
    if (!convert_integer_data(data, integer_data, integer_error)) {
        return fail_model(integer_error);
    }
    std::int64_t route_horizon = integer_data.limit;
    if (options.solve_mode == CpSolveMode::ThresholdOptimization &&
        !compute_threshold_horizon(
            integer_data.limit,
            options.threshold_horizon_factor,
            route_horizon,
            integer_error
        )) {
        return fail_model(integer_error);
    }
    result.model_horizon = route_horizon;

    const int n = static_cast<int>(data.requests.size());
    const int k_count = options.vehicle_count;
    if (n == 0 || k_count > n) {
        result.outcome = CpSolveOutcome::ProvenInfeasible;
        result.status_name = "INFEASIBLE";
        result.termination_name = "proven-infeasible";
        return result;
    }
    if (graph.end_node_id() != 2 * n + 1) {
        return fail_model("Multigraph service-node indexing is inconsistent with the requests.");
    }
    if (route_horizon >
        std::numeric_limits<std::int64_t>::max() /
            std::max<std::int64_t>(2, k_count)) {
        return fail_model("Active CP horizon makes a redundant bound overflow.");
    }

    const NodeId end_node = graph.end_node_id();
    const int service_count = 2 * n;
    auto service_index = [service_count](NodeId node) {
        return node >= 1 && node <= service_count ? node - 1 : -1;
    };
    auto start_anchor = [service_count](int k) {
        return service_count + 2 * k;
    };
    auto end_anchor = [service_count](int k) {
        return service_count + 2 * k + 1;
    };

    for (NodeId node = 1; node <= service_count; ++node) {
        try {
            const NodeSpec& spec = graph.node(node);
            if (spec.kind != NodeSpec::Kind::Pickup &&
                spec.kind != NodeSpec::Kind::Delivery) {
                return fail_model("A physical multigraph node is not pickup or delivery.");
            }
            if (!spec.request_idx.has_value() ||
                *spec.request_idx < 0 || *spec.request_idx >= n) {
                return fail_model("A physical multigraph node has no valid request index.");
            }
        } catch (const std::exception& error) {
            return fail_model(error.what());
        }
    }

    std::vector<std::int64_t> edge_time(graph.number_of_edges(), 0);
    std::vector<int> start_edge_ids;
    std::vector<int> internal_edge_ids;
    std::vector<int> end_edge_ids;
    for (std::size_t edge_id = 0; edge_id < graph.edges().size(); ++edge_id) {
        const EdgeRecord& edge = graph.edges()[edge_id];
        if (!exact_integer(edge.data.time, edge_time[edge_id]) ||
            edge_time[edge_id] < 0) {
            return fail_model(
                "Multigraph CP-SAT requires nonnegative exact integer edge durations."
            );
        }
        if (edge.u == 0 && edge.v == end_node) {
            continue;
        }
        if (edge.u == 0 && service_index(edge.v) >= 0 &&
            graph.node(edge.v).kind == NodeSpec::Kind::Pickup) {
            start_edge_ids.push_back(static_cast<int>(edge_id));
        } else if (service_index(edge.u) >= 0 &&
                   service_index(edge.v) >= 0) {
            internal_edge_ids.push_back(static_cast<int>(edge_id));
        } else if (service_index(edge.u) >= 0 && edge.v == end_node &&
                   graph.node(edge.u).kind == NodeSpec::Kind::Delivery) {
            end_edge_ids.push_back(static_cast<int>(edge_id));
        } else {
            return fail_model("Multigraph contains an unsupported non-dummy edge.");
        }
    }
    if (start_edge_ids.empty() || end_edge_ids.empty()) {
        result.outcome = CpSolveOutcome::ProvenInfeasible;
        result.status_name = "INFEASIBLE";
        result.termination_name = "proven-infeasible";
        return result;
    }

    CpModelBuilder builder;
    CircuitConstraint circuit = builder.AddCircuitConstraint();
    std::vector<MultigraphArcVariable> selected_arcs;
    selected_arcs.reserve(
        internal_edge_ids.size() +
        static_cast<std::size_t>(k_count) *
            (start_edge_ids.size() + end_edge_ids.size())
    );
    std::vector<std::vector<MultigraphArcVariable>> start_arcs_by_vehicle(
        static_cast<std::size_t>(k_count)
    );
    std::vector<std::map<State, StateFlowTerms>> state_flow(
        static_cast<std::size_t>(end_node + 1)
    );

    for (int edge_id : internal_edge_ids) {
        const EdgeRecord& edge = graph.edges()[static_cast<std::size_t>(edge_id)];
        BoolVar literal = builder.NewBoolVar().WithName(
            "mg_arc_e" + std::to_string(edge_id)
        );
        circuit.AddArc(service_index(edge.u), service_index(edge.v), literal);
        MultigraphArcVariable arc{
            service_index(edge.u), service_index(edge.v), edge_id, -1, literal
        };
        selected_arcs.push_back(arc);
        state_flow[static_cast<std::size_t>(edge.u)][edge.data.start_state]
            .outgoing.push_back(literal);
        state_flow[static_cast<std::size_t>(edge.v)][edge.data.end_state]
            .incoming.push_back(literal);
        ++result.build_stats.created_multigraph_internal_arcs;
    }

    for (int k = 0; k < k_count; ++k) {
        for (int edge_id : start_edge_ids) {
            const EdgeRecord& edge = graph.edges()[static_cast<std::size_t>(edge_id)];
            BoolVar literal = builder.NewBoolVar().WithName(
                "mg_arc_S" + std::to_string(k) + "_e" + std::to_string(edge_id)
            );
            circuit.AddArc(start_anchor(k), service_index(edge.v), literal);
            MultigraphArcVariable arc{
                start_anchor(k), service_index(edge.v), edge_id, k, literal
            };
            selected_arcs.push_back(arc);
            start_arcs_by_vehicle[static_cast<std::size_t>(k)].push_back(arc);
            state_flow[static_cast<std::size_t>(edge.v)][edge.data.end_state]
                .incoming.push_back(literal);
            ++result.build_stats.created_multigraph_start_arc_copies;
        }
        for (int edge_id : end_edge_ids) {
            const EdgeRecord& edge = graph.edges()[static_cast<std::size_t>(edge_id)];
            BoolVar literal = builder.NewBoolVar().WithName(
                "mg_arc_e" + std::to_string(edge_id) + "_Z" + std::to_string(k)
            );
            circuit.AddArc(service_index(edge.u), end_anchor(k), literal);
            selected_arcs.push_back(MultigraphArcVariable{
                service_index(edge.u), end_anchor(k), edge_id, k, literal
            });
            state_flow[static_cast<std::size_t>(edge.u)][edge.data.start_state]
                .outgoing.push_back(literal);
            ++result.build_stats.created_multigraph_end_arc_copies;
        }
        circuit.AddArc(
            end_anchor(k),
            start_anchor((k + 1) % k_count),
            builder.TrueVar()
        );
        ++result.build_stats.fixed_multigraph_connector_arcs;
    }

    for (NodeId node = 1; node <= service_count; ++node) {
        for (const auto& entry : state_flow[static_cast<std::size_t>(node)]) {
            builder.AddEquality(
                bool_sum(entry.second.incoming),
                bool_sum(entry.second.outgoing)
            );
            ++result.build_stats.state_continuity_constraint_count;
        }
    }

    std::vector<IntVar> vehicle;
    std::vector<IntVar> position;
    std::vector<IntVar> after_service;
    std::vector<std::vector<BoolVar>> assigned(
        static_cast<std::size_t>(service_count),
        std::vector<BoolVar>(static_cast<std::size_t>(k_count))
    );
    vehicle.reserve(static_cast<std::size_t>(service_count));
    position.reserve(static_cast<std::size_t>(service_count));
    after_service.reserve(static_cast<std::size_t>(service_count));
    for (int index = 0; index < service_count; ++index) {
        vehicle.push_back(builder.NewIntVar(Domain(0, k_count - 1)).WithName(
            "mg_vehicle_" + std::to_string(index + 1)
        ));
        position.push_back(builder.NewIntVar(Domain(1, service_count)).WithName(
            "mg_position_" + std::to_string(index + 1)
        ));
        after_service.push_back(builder.NewIntVar(Domain(0, route_horizon)).WithName(
            "mg_after_" + std::to_string(index + 1)
        ));
        std::vector<BoolVar> assignment_row;
        for (int k = 0; k < k_count; ++k) {
            BoolVar literal = builder.NewBoolVar().WithName(
                "mg_assigned_" + std::to_string(index + 1) + "_" +
                std::to_string(k)
            );
            assigned[static_cast<std::size_t>(index)][static_cast<std::size_t>(k)] =
                literal;
            builder.AddEquality(vehicle[static_cast<std::size_t>(index)], k)
                .OnlyEnforceIf(literal);
            assignment_row.push_back(literal);
        }
        builder.AddExactlyOne(assignment_row);
    }

    for (const MultigraphArcVariable& arc : selected_arcs) {
        const EdgeRecord& edge =
            graph.edges()[static_cast<std::size_t>(arc.original_edge_id)];
        const std::int64_t duration =
            edge_time[static_cast<std::size_t>(arc.original_edge_id)];
        if (edge.u == 0) {
            const int head = service_index(edge.v);
            builder.AddImplication(
                arc.literal,
                assigned[static_cast<std::size_t>(head)]
                        [static_cast<std::size_t>(arc.vehicle_index)]
            );
            builder.AddEquality(position[static_cast<std::size_t>(head)], 1)
                .OnlyEnforceIf(arc.literal);
            builder.AddEquality(after_service[static_cast<std::size_t>(head)], duration)
                .OnlyEnforceIf(arc.literal);
        } else if (edge.v == end_node) {
            const int tail = service_index(edge.u);
            builder.AddImplication(
                arc.literal,
                assigned[static_cast<std::size_t>(tail)]
                        [static_cast<std::size_t>(arc.vehicle_index)]
            );
        } else {
            const int tail = service_index(edge.u);
            const int head = service_index(edge.v);
            builder.AddEquality(
                vehicle[static_cast<std::size_t>(head)],
                vehicle[static_cast<std::size_t>(tail)]
            ).OnlyEnforceIf(arc.literal);
            builder.AddEquality(
                position[static_cast<std::size_t>(head)],
                position[static_cast<std::size_t>(tail)] + 1
            ).OnlyEnforceIf(arc.literal);
            builder.AddEquality(
                after_service[static_cast<std::size_t>(head)],
                after_service[static_cast<std::size_t>(tail)] + duration
            ).OnlyEnforceIf(arc.literal);
        }
    }

    std::vector<IntVar> completion;
    completion.reserve(static_cast<std::size_t>(k_count));
    for (int k = 0; k < k_count; ++k) {
        completion.push_back(builder.NewIntVar(Domain(0, route_horizon)).WithName(
            "mg_completion_" + std::to_string(k)
        ));
    }
    for (const MultigraphArcVariable& arc : selected_arcs) {
        const EdgeRecord& edge =
            graph.edges()[static_cast<std::size_t>(arc.original_edge_id)];
        if (edge.v != end_node) {
            continue;
        }
        builder.AddEquality(
            completion[static_cast<std::size_t>(arc.vehicle_index)],
            after_service[static_cast<std::size_t>(service_index(edge.u))] +
                edge_time[static_cast<std::size_t>(arc.original_edge_id)]
        ).OnlyEnforceIf(arc.literal);
    }

    std::optional<IntVar> makespan;
    if (options.solve_mode == CpSolveMode::ThresholdOptimization) {
        makespan.emplace(
            builder.NewIntVar(Domain(0, route_horizon)).WithName("mg_makespan")
        );
        builder.AddMaxEquality(*makespan, completion);
        builder.Minimize(*makespan);
    }

    std::set<int> types;
    for (const Request& request : data.requests) {
        types.insert(request.container_type);
    }
    if (options.redundant.terminal_balance) {
        for (int k = 0; k < k_count; ++k) {
            LinearExpr pickups;
            LinearExpr deliveries;
            for (int i = 0; i < n; ++i) {
                pickups += assigned[static_cast<std::size_t>(i)]
                                   [static_cast<std::size_t>(k)];
                deliveries += assigned[static_cast<std::size_t>(n + i)]
                                      [static_cast<std::size_t>(k)];
            }
            builder.AddEquality(pickups, deliveries);
            ++result.build_stats.terminal_balance_constraint_count;
            for (int type : types) {
                LinearExpr typed_pickups;
                LinearExpr typed_deliveries;
                for (int i = 0; i < n; ++i) {
                    if (data.requests[static_cast<std::size_t>(i)].container_type != type) {
                        continue;
                    }
                    typed_pickups += assigned[static_cast<std::size_t>(i)]
                                             [static_cast<std::size_t>(k)];
                    typed_deliveries += assigned[static_cast<std::size_t>(n + i)]
                                                [static_cast<std::size_t>(k)];
                }
                builder.AddEquality(typed_pickups, typed_deliveries);
                ++result.build_stats.terminal_balance_constraint_count;
            }
        }
    }

    if (options.redundant.container_workload) {
        for (int k = 0; k < k_count; ++k) {
            LinearExpr service_twice;
            LinearExpr workload;
            std::map<std::pair<int, int>, BoolVar> matching;
            for (int i = 0; i < n; ++i) {
                const Request& request = data.requests[static_cast<std::size_t>(i)];
                service_twice += 2 * (integer_data.pickup + integer_data.treatment) *
                    assigned[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
                service_twice += 2 * integer_data.delivery *
                    assigned[static_cast<std::size_t>(n + i)][static_cast<std::size_t>(k)];
                workload += integer_data.travel[static_cast<std::size_t>(request.from_id)]
                                               [static_cast<std::size_t>(request.to_id)] *
                    assigned[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
                std::vector<BoolVar> outgoing;
                for (int j = 0; j < n; ++j) {
                    if (request.container_type !=
                        data.requests[static_cast<std::size_t>(j)].container_type) {
                        continue;
                    }
                    BoolVar lambda = builder.NewBoolVar().WithName(
                        "mg_lambda_" + std::to_string(i) + "_" +
                        std::to_string(j) + "_" + std::to_string(k)
                    );
                    matching.emplace(std::make_pair(i, j), lambda);
                    outgoing.push_back(lambda);
                    workload += integer_data.travel[static_cast<std::size_t>(request.to_id)]
                                                   [static_cast<std::size_t>(
                                                       data.requests[static_cast<std::size_t>(j)]
                                                           .from_id)] * lambda;
                }
                builder.AddEquality(
                    bool_sum(outgoing),
                    assigned[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]
                );
                ++result.build_stats.container_workload_constraint_count;
            }
            for (int j = 0; j < n; ++j) {
                std::vector<BoolVar> incoming;
                for (int i = 0; i < n; ++i) {
                    const auto found = matching.find(std::make_pair(i, j));
                    if (found != matching.end()) {
                        incoming.push_back(found->second);
                    }
                }
                builder.AddEquality(
                    bool_sum(incoming),
                    assigned[static_cast<std::size_t>(n + j)]
                            [static_cast<std::size_t>(k)]
                );
                ++result.build_stats.container_workload_constraint_count;
            }
            builder.AddLessOrEqual(
                service_twice + workload,
                2 * route_horizon
            );
            ++result.build_stats.container_workload_constraint_count;
        }
    }

    if (options.redundant.aggregate_duration) {
        LinearExpr total_duration;
        for (const MultigraphArcVariable& arc : selected_arcs) {
            total_duration +=
                edge_time[static_cast<std::size_t>(arc.original_edge_id)] * arc.literal;
        }
        builder.AddLessOrEqual(total_duration, k_count * route_horizon);
        result.build_stats.aggregate_duration_constraint_count = 1;
    }

    if (options.symmetry.first_pickup_vehicle_ordering) {
        for (int k = 0; k + 1 < k_count; ++k) {
            LinearExpr current;
            LinearExpr next;
            for (const MultigraphArcVariable& arc :
                 start_arcs_by_vehicle[static_cast<std::size_t>(k)]) {
                const EdgeRecord& edge =
                    graph.edges()[static_cast<std::size_t>(arc.original_edge_id)];
                current += (*graph.node(edge.v).request_idx + 1) * arc.literal;
            }
            for (const MultigraphArcVariable& arc :
                 start_arcs_by_vehicle[static_cast<std::size_t>(k + 1)]) {
                const EdgeRecord& edge =
                    graph.edges()[static_cast<std::size_t>(arc.original_edge_id)];
                next += (*graph.node(edge.v).request_idx + 1) * arc.literal;
            }
            builder.AddLessOrEqual(current + 1, next);
            ++result.build_stats.first_pickup_symmetry_constraint_count;
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
                        after_service[static_cast<std::size_t>(requests[left])],
                        after_service[static_cast<std::size_t>(requests[right])]
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
    parameters.set_log_search_progress(write_log_file);
    parameters.set_log_to_stdout(false);
    operations_research::sat::Model model;
    model.Add(operations_research::sat::NewSatParameters(parameters));
    if (write_log_file) {
        operations_research::SolverLogger* logger =
            model.GetOrCreate<operations_research::SolverLogger>();
        logger->EnableLogging(true);
        logger->SetLogToStdOut(false);
        logger->AddInfoLoggingCallback([&log_file](const std::string& message) {
            log_file << message;
            if (message.empty() || message.back() != '\n') {
                log_file << '\n';
            }
            log_file.flush();
        });
    }

    std::mutex callback_mutex;
    std::optional<CpSolverResponse> accepted_threshold_response;
    bool stopped_by_bound_callback = false;
    double callback_proving_bound = 0.0;
    if (options.solve_mode == CpSolveMode::ThresholdOptimization) {
        model.Add(operations_research::sat::NewFeasibleSolutionObserver(
            [&](const CpSolverResponse& candidate) {
                if (candidate.objective_value() >
                    static_cast<double>(integer_data.limit) + kObjectiveTolerance) {
                    return;
                }
                bool should_stop = false;
                {
                    std::lock_guard<std::mutex> lock(callback_mutex);
                    if (!accepted_threshold_response.has_value()) {
                        accepted_threshold_response = candidate;
                        should_stop = true;
                    }
                }
                if (should_stop) {
                    operations_research::sat::StopSearch(&model);
                }
            }
        ));
        model.Add(operations_research::sat::NewBestBoundCallback(
            [&](double bound) {
                if (!objective_bound_excludes_threshold(bound, integer_data.limit)) {
                    return;
                }
                bool should_stop = false;
                {
                    std::lock_guard<std::mutex> lock(callback_mutex);
                    callback_proving_bound =
                        std::max(callback_proving_bound, bound);
                    if (!stopped_by_bound_callback) {
                        stopped_by_bound_callback = true;
                        should_stop = true;
                    }
                }
                if (should_stop) {
                    operations_research::sat::StopSearch(&model);
                }
            }
        ));
    }

    const CpSolverResponse response =
        operations_research::sat::SolveCpModel(builder.Build(), &model);
    std::optional<CpSolverResponse> saved_threshold_response;
    {
        std::lock_guard<std::mutex> lock(callback_mutex);
        saved_threshold_response = accepted_threshold_response;
        result.stopped_by_feasible_observer =
            accepted_threshold_response.has_value();
        result.stopped_by_bound_callback = stopped_by_bound_callback;
        if (stopped_by_bound_callback) {
            result.has_objective_bound = true;
            result.best_objective_bound = callback_proving_bound;
        }
    }

    result.raw_status = static_cast<int>(response.status());
    result.status_name =
        operations_research::sat::CpSolverStatus_Name(response.status());
    result.solution_info = response.solution_info();
    result.wall_time_seconds = response.wall_time();
    result.conflicts = response.num_conflicts();
    result.branches = response.num_branches();
    const bool raw_has_solution =
        response.status() == operations_research::sat::CpSolverStatus::FEASIBLE ||
        response.status() == operations_research::sat::CpSolverStatus::OPTIMAL;
    if (options.solve_mode == CpSolveMode::ThresholdOptimization) {
        if (raw_has_solution) {
            result.has_objective_value = true;
            result.objective_value = response.objective_value();
        }
        if (std::isfinite(response.best_objective_bound())) {
            result.has_objective_bound = true;
            result.best_objective_bound = std::max(
                result.best_objective_bound,
                response.best_objective_bound()
            );
        }
        if (saved_threshold_response.has_value()) {
            result.has_objective_value = true;
            result.objective_value = saved_threshold_response->objective_value();
        }
    }

    auto append_final_log = [&]() {
        if (!write_log_file) {
            return;
        }
        log_file << "[spdp-cp-wrapper] graph_mode="
                 << cp_graph_mode_name(result.graph_mode)
                 << " mode=" << cp_solve_mode_name(result.solve_mode)
                 << " horizon_factor=" << result.threshold_horizon_factor
                 << " model_horizon=" << result.model_horizon
                 << " termination=" << result.termination_name
                 << " objective_value=" << result.objective_value
                 << " best_objective_bound=" << result.best_objective_bound
                 << " full_skip_reservoir_requested="
                 << (options.redundant.full_skip_reservoir ? 1 : 0)
                 << " full_skip_reservoir_effective=embedded-state"
                 << " internal_arcs="
                 << result.build_stats.created_multigraph_internal_arcs
                 << " start_arc_copies="
                 << result.build_stats.created_multigraph_start_arc_copies
                 << " end_arc_copies="
                 << result.build_stats.created_multigraph_end_arc_copies
                 << '\n';
        log_file.flush();
    };

    if (response.status() == operations_research::sat::CpSolverStatus::MODEL_INVALID) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.termination_name = "model-invalid";
        result.error_message = result.solution_info;
        append_final_log();
        return result;
    }

    const CpSolverResponse* solution_response = nullptr;
    if (options.solve_mode == CpSolveMode::Satisfaction) {
        if (response.status() == operations_research::sat::CpSolverStatus::INFEASIBLE) {
            result.outcome = CpSolveOutcome::ProvenInfeasible;
            result.termination_name = "proven-infeasible";
            append_final_log();
            return result;
        }
        if (raw_has_solution) {
            result.outcome = CpSolveOutcome::Feasible;
            result.termination_name = "feasible";
            solution_response = &response;
        }
    } else {
        if (saved_threshold_response.has_value()) {
            result.outcome = CpSolveOutcome::Feasible;
            result.termination_name = "threshold-feasible-witness";
            solution_response = &*saved_threshold_response;
        } else if (raw_has_solution &&
                   response.objective_value() <=
                       static_cast<double>(integer_data.limit) + kObjectiveTolerance) {
            result.outcome = CpSolveOutcome::Feasible;
            result.termination_name = "threshold-feasible-witness";
            solution_response = &response;
        } else if (result.stopped_by_bound_callback ||
                   (result.has_objective_bound && objective_bound_excludes_threshold(
                       result.best_objective_bound, integer_data.limit)) ||
                   (response.status() ==
                        operations_research::sat::CpSolverStatus::OPTIMAL &&
                    raw_has_solution &&
                    response.objective_value() >
                        static_cast<double>(integer_data.limit) + kObjectiveTolerance)) {
            result.outcome = CpSolveOutcome::ProvenInfeasible;
            result.termination_name = "objective-bound-infeasible";
            append_final_log();
            return result;
        } else if (response.status() ==
                   operations_research::sat::CpSolverStatus::INFEASIBLE) {
            result.outcome = CpSolveOutcome::ProvenInfeasible;
            result.termination_name = "restricted-horizon-infeasible";
            append_final_log();
            return result;
        }
    }

    if (solution_response == nullptr) {
        result.outcome = classify_cp_unknown_outcome(
            result.configured_time_limit_seconds,
            result.wall_time_seconds
        );
        result.hit_time_limit =
            result.outcome == CpSolveOutcome::TimedOutUnknown;
        result.early_unknown = result.outcome == CpSolveOutcome::EarlyUnknown;
        result.termination_name = result.hit_time_limit
            ? "timed-out-unknown" : "early-unknown";
        append_final_log();
        return result;
    }

    std::set<int> seen_edge_ids;
    for (const MultigraphArcVariable& arc : selected_arcs) {
        if (!operations_research::sat::SolutionBooleanValue(
                *solution_response, arc.literal)) {
            continue;
        }
        if (!seen_edge_ids.insert(arc.original_edge_id).second) {
            result.outcome = CpSolveOutcome::ModelInvalid;
            result.termination_name = "model-invalid";
            result.error_message =
                "Multigraph CP witness selects one original edge more than once.";
            result.active_edge_ids.clear();
            append_final_log();
            return result;
        }
        result.active_edge_ids.push_back(arc.original_edge_id);
    }
    if (result.active_edge_ids.size() !=
        static_cast<std::size_t>(service_count + k_count)) {
        result.outcome = CpSolveOutcome::ModelInvalid;
        result.termination_name = "model-invalid";
        result.error_message =
            "Multigraph CP witness has an unexpected selected-edge count.";
        result.active_edge_ids.clear();
    }
    append_final_log();
    return result;
}

CpFixedKSolveResult solve_fixed_k_cp_sat(
    const SPDPData& data,
    const CpSatSolveOptions& options
) {
    CpSatSolveOptions original_options = options;
    original_options.graph_mode = CpGraphMode::OriginalGraph;
    return solve_fixed_k_original_graph_cp_sat(data, original_options);
}

CpFixedKSolveResult solve_fixed_k_cp_sat(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const CpSatSolveOptions& options
) {
    switch (options.graph_mode) {
    case CpGraphMode::OriginalGraph:
        return solve_fixed_k_original_graph_cp_sat(data, options);
    case CpGraphMode::Multigraph:
        return solve_fixed_k_multigraph_cp_sat(data, graph, options);
    }
    CpFixedKSolveResult result;
    result.outcome = CpSolveOutcome::ModelInvalid;
    result.status_name = "MODEL_INVALID";
    result.termination_name = "model-invalid";
    result.error_message = "Unsupported CP graph mode.";
    return result;
}

}  // namespace spdp
