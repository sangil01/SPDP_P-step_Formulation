#include "KMinComputation.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "TwoIndexModelCore.h"
#include "gurobi_c++.h"

namespace spdp {
namespace {

constexpr double kTolerance = 1e-9;
constexpr double kAuxiliaryLPSolverTolerance = 1e-9;

double safe_duration_lower_bound(
    double objective_bound,
    double route_time_limit
) {
    const double scale = std::max(
        {1.0, route_time_limit, std::abs(objective_bound)}
    );
    const double tolerance = 10.0 * kAuxiliaryLPSolverTolerance * scale;
    return std::max(0.0, objective_bound - tolerance);
}

using detail::TwoIndexCoreModel;
using detail::TwoIndexCoreObjective;
using detail::TwoIndexCoreOptions;
using detail::build_two_index_model_core;

int cor_route_lower_bound_kmin(const SPDPData& data) {
    if (data.time_limit <= 0.0) {
        throw std::runtime_error("Time limit must be positive for VI-44.");
    }

    double pickup_to_treatment_time_sum = 0.0;
    double treatment_to_best_delivery_time_sum = 0.0;
    for (const Request& request : data.requests) {
        pickup_to_treatment_time_sum +=
            data.time[static_cast<std::size_t>(request.from_id)]
                     [static_cast<std::size_t>(request.to_id)];

        double best_compatible_delivery_time = std::numeric_limits<double>::infinity();
        for (const Request& candidate : data.requests) {
            if (candidate.container_type != request.container_type) {
                continue;
            }
            best_compatible_delivery_time = std::min(
                best_compatible_delivery_time,
                data.time[static_cast<std::size_t>(request.to_id)]
                         [static_cast<std::size_t>(candidate.from_id)]
            );
        }
        if (!std::isfinite(best_compatible_delivery_time)) {
            throw std::runtime_error(
                "Failed to compute VI-44 k_min because no compatible delivery location exists."
            );
        }
        treatment_to_best_delivery_time_sum += best_compatible_delivery_time;
    }

    const double total_service_time =
        static_cast<double>(data.requests.size()) *
        (data.time_pickup + data.time_empty + data.time_delivery);
    const double lower_bound =
        0.5 * (pickup_to_treatment_time_sum + treatment_to_best_delivery_time_sum) +
        total_service_time;
    return std::max(
        0,
        static_cast<int>(std::ceil(lower_bound / data.time_limit - kTolerance))
    );
}

struct AuxiliaryDurationSubproblemResult {
    int status = 0;
    bool hit_time_limit = false;
    bool has_certified_bound = false;
    double objective_value = -1.0;
    double objective_bound = -1.0;
    double runtime_seconds = 0.0;
    bool stopped_by_rounded_bound = false;
    bool rounded_bound_certified = false;
    int certified_rounded_k = 0;
    double callback_objective_ub = -1.0;
    double callback_safe_objective_lb = -1.0;
};

class RoundedDurationBoundCallback final : public GRBCallback {
public:
    RoundedDurationBoundCallback(
        const MultiDiGraph& graph,
        const std::vector<GRBVar>& y_vars,
        double route_time_limit
    )
        : graph_(graph),
          y_vars_(y_vars),
          route_time_limit_(route_time_limit) {}

    bool stopped() const { return stopped_; }
    int certified_k() const { return certified_k_; }
    double exact_upper_bound() const { return exact_upper_bound_; }
    double safe_lower_bound() const { return safe_lower_bound_; }

protected:
    void callback() override {
        try {
            if (where == GRB_CB_MIPSOL) {
                double exact_objective = 0.0;
                for (std::size_t edge_id = 0; edge_id < y_vars_.size(); ++edge_id) {
                    if (getSolution(y_vars_[edge_id]) > 0.5) {
                        exact_objective += graph_.edges()[edge_id].data.time;
                    }
                }
                exact_upper_bound_ = exact_objective;
                maybe_stop(getDoubleInfo(GRB_CB_MIPSOL_OBJBND));
                return;
            }
            if (where == GRB_CB_MIP && exact_upper_bound_ >= 0.0) {
                maybe_stop(getDoubleInfo(GRB_CB_MIP_OBJBND));
            }
        } catch (const GRBException&) {
            // Continue safely; post-solve attributes still give a valid result.
        }
    }

private:
    void maybe_stop(double objective_bound) {
        const std::optional<int> certificate = certified_rounded_duration_bound(
            exact_upper_bound_, objective_bound, route_time_limit_);
        if (!certificate.has_value()) {
            return;
        }
        stopped_ = true;
        certified_k_ = certificate.value();
        safe_lower_bound_ = safe_duration_lower_bound(
            objective_bound, route_time_limit_);
        abort();
    }

    const MultiDiGraph& graph_;
    const std::vector<GRBVar>& y_vars_;
    double route_time_limit_ = 0.0;
    bool stopped_ = false;
    int certified_k_ = 0;
    double exact_upper_bound_ = -1.0;
    double safe_lower_bound_ = -1.0;
};

AuxiliaryDurationSubproblemResult solve_auxiliary_duration_subproblem(
    const SPDPData& data,
    const MultiDiGraph& graph,
    double solver_time_limit,
    bool binary_y,
    bool add_time_constraints,
    bool rounded_bound_stop
) {
    TwoIndexCoreOptions core_options;
    core_options.objective = TwoIndexCoreObjective::Duration;
    core_options.binary_y = binary_y;
    core_options.add_time_constraints = add_time_constraints;
    core_options.solver_time_limit = solver_time_limit;
    core_options.gurobi_threads = 1;
    core_options.output_enabled = false;
    core_options.name_prefix = "vi44_aux";
    TwoIndexCoreModel core = build_two_index_model_core(data, graph, core_options);

    std::unique_ptr<RoundedDurationBoundCallback> callback;
    if (binary_y && rounded_bound_stop) {
        callback = std::make_unique<RoundedDurationBoundCallback>(
            graph, core.y_vars, data.time_limit);
        core.model->setCallback(callback.get());
    }

    const auto solve_start = std::chrono::steady_clock::now();
    core.model->optimize();
    const auto solve_end = std::chrono::steady_clock::now();

    AuxiliaryDurationSubproblemResult result;
    result.status = core.model->get(GRB_IntAttr_Status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.runtime_seconds =
        std::chrono::duration<double>(solve_end - solve_start).count();

    if (core.model->get(GRB_IntAttr_SolCount) > 0) {
        try {
            const double objective_value = core.model->get(GRB_DoubleAttr_ObjVal);
            if (std::isfinite(objective_value) &&
                std::abs(objective_value) < 0.5 * GRB_INFINITY) {
                result.objective_value = objective_value;
            }
        } catch (const GRBException& error) {
            if (error.getErrorCode() != GRB_ERROR_DATA_NOT_AVAILABLE) {
                throw;
            }
        }
    }

    try {
        const double objective_bound = core.model->get(GRB_DoubleAttr_ObjBound);
        if (std::isfinite(objective_bound) &&
            std::abs(objective_bound) < 0.5 * GRB_INFINITY) {
            result.objective_bound = objective_bound;
            result.has_certified_bound = true;
        }
    } catch (const GRBException& error) {
        if (error.getErrorCode() != GRB_ERROR_DATA_NOT_AVAILABLE) {
            throw;
        }
    }

    if (callback != nullptr && callback->stopped()) {
        result.stopped_by_rounded_bound = true;
        result.certified_rounded_k = callback->certified_k();
        result.callback_objective_ub = callback->exact_upper_bound();
        result.callback_safe_objective_lb = callback->safe_lower_bound();
        const std::optional<int> postsolve_certificate =
            result.has_certified_bound
                ? certified_rounded_duration_bound(
                      result.callback_objective_ub,
                      result.objective_bound,
                      data.time_limit)
                : std::nullopt;
        result.rounded_bound_certified =
            postsolve_certificate.has_value() &&
            postsolve_certificate.value() == result.certified_rounded_k;
    }

    const bool supported_interrupted =
        result.status == GRB_INTERRUPTED &&
        result.stopped_by_rounded_bound &&
        result.rounded_bound_certified;
    if (result.status != GRB_OPTIMAL && result.status != GRB_TIME_LIMIT &&
        !supported_interrupted) {
        throw std::runtime_error(
            "The auxiliary VI-44 subproblem stopped with an unsupported status (" +
            std::to_string(result.status) + ")."
        );
    }
    if (result.status == GRB_OPTIMAL && !result.has_certified_bound) {
        throw std::runtime_error(
            "The optimal auxiliary VI-44 subproblem returned no finite certified bound."
        );
    }
    return result;
}

struct VehicleAssignmentMatchingArc {
    std::size_t pickup_index = 0U;
    std::size_t delivery_index = 0U;
    std::vector<GRBVar> variables;
};

struct VehicleAssignmentPhysicalArc {
    int tail = -1;
    int head = -1;
    GRBVar selected;
    GRBVar flow;
};

struct VehicleAssignmentFeasibilityResult {
    int status = 0;
    bool hit_time_limit = false;
    bool feasible = false;
    bool proven_infeasible = false;
    double runtime_seconds = 0.0;
    std::vector<VI44VehicleAssignmentVehicleResult> vehicles;
};

VehicleAssignmentFeasibilityResult solve_vehicle_assignment_feasibility(
    const SPDPData& data,
    int vehicle_count,
    const VI44VehicleAssignmentOptions& options,
    GRBEnv& environment,
    double solver_time_limit
) {
    const std::size_t request_count = data.requests.size();
    const int location_count = data.locations;
    if (vehicle_count <= 0 ||
        static_cast<std::size_t>(vehicle_count) > request_count) {
        throw std::runtime_error(
            "The VI-44 vehicle-assignment vehicle count must be in [1, |I|]."
        );
    }
    if (location_count <= 0) {
        throw std::runtime_error(
            "The VI-44 vehicle-assignment problem requires physical locations."
        );
    }

    GRBModel model(environment);
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_Threads, 1);
    model.set(GRB_IntParam_DualReductions, 0);
    model.set(GRB_IntParam_SolutionLimit, 1);
    model.set(GRB_DoubleParam_FeasibilityTol, kAuxiliaryLPSolverTolerance);
    model.set(GRB_DoubleParam_OptimalityTol, kAuxiliaryLPSolverTolerance);
    if (solver_time_limit > 0.0) {
        model.set(GRB_DoubleParam_TimeLimit, solver_time_limit);
    }

    std::vector<std::vector<GRBVar>> pickup_assignment(request_count);
    std::vector<std::vector<GRBVar>> delivery_assignment(request_count);
    for (std::size_t request_index = 0U; request_index < request_count; ++request_index) {
        pickup_assignment[request_index].reserve(static_cast<std::size_t>(vehicle_count));
        delivery_assignment[request_index].reserve(static_cast<std::size_t>(vehicle_count));
        for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
            pickup_assignment[request_index].push_back(model.addVar(
                0.0,
                1.0,
                0.0,
                GRB_BINARY,
                "vi44_va_u_" + std::to_string(request_index) + "_" +
                    std::to_string(vehicle)
            ));
            delivery_assignment[request_index].push_back(model.addVar(
                0.0,
                1.0,
                0.0,
                GRB_BINARY,
                "vi44_va_v_" + std::to_string(request_index) + "_" +
                    std::to_string(vehicle)
            ));
        }
    }

    std::vector<VehicleAssignmentMatchingArc> matching_arcs;
    for (std::size_t pickup_index = 0U;
         pickup_index < request_count;
         ++pickup_index) {
        for (std::size_t delivery_index = 0U;
             delivery_index < request_count;
             ++delivery_index) {
            if (data.requests[pickup_index].container_type !=
                data.requests[delivery_index].container_type) {
                continue;
            }
            VehicleAssignmentMatchingArc arc;
            arc.pickup_index = pickup_index;
            arc.delivery_index = delivery_index;
            arc.variables.reserve(static_cast<std::size_t>(vehicle_count));
            for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
                arc.variables.push_back(model.addVar(
                    0.0,
                    1.0,
                    0.0,
                    GRB_BINARY,
                    "vi44_va_mu_" + std::to_string(pickup_index) + "_" +
                        std::to_string(delivery_index) + "_" +
                        std::to_string(vehicle)
                ));
            }
            matching_arcs.push_back(std::move(arc));
        }
    }

    std::vector<std::vector<GRBVar>> required_location;
    std::vector<std::vector<GRBVar>> path_start;
    std::vector<std::vector<GRBVar>> path_end;
    std::vector<std::vector<GRBVar>> source_flow;
    std::vector<std::vector<VehicleAssignmentPhysicalArc>> physical_arcs;
    if (options.add_tsp_bound) {
        required_location.resize(static_cast<std::size_t>(location_count));
        path_start.resize(static_cast<std::size_t>(location_count));
        path_end.resize(static_cast<std::size_t>(location_count));
        source_flow.resize(static_cast<std::size_t>(location_count));
        for (int location = 0; location < location_count; ++location) {
            for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
                required_location[static_cast<std::size_t>(location)].push_back(
                    model.addVar(
                        0.0,
                        1.0,
                        0.0,
                        GRB_BINARY,
                        "vi44_va_b_" + std::to_string(location) + "_" +
                            std::to_string(vehicle)
                    )
                );
                path_start[static_cast<std::size_t>(location)].push_back(
                    model.addVar(
                        0.0,
                        1.0,
                        0.0,
                        GRB_BINARY,
                        "vi44_va_rho_" + std::to_string(location) + "_" +
                            std::to_string(vehicle)
                    )
                );
                path_end[static_cast<std::size_t>(location)].push_back(
                    model.addVar(
                        0.0,
                        1.0,
                        0.0,
                        GRB_BINARY,
                        "vi44_va_eta_" + std::to_string(location) + "_" +
                            std::to_string(vehicle)
                    )
                );
                source_flow[static_cast<std::size_t>(location)].push_back(
                    model.addVar(
                        0.0,
                        static_cast<double>(location_count),
                        0.0,
                        GRB_CONTINUOUS,
                        "vi44_va_fs_" + std::to_string(location) + "_" +
                            std::to_string(vehicle)
                    )
                );
            }
        }

        physical_arcs.resize(static_cast<std::size_t>(vehicle_count));
        for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
            auto& vehicle_arcs = physical_arcs[static_cast<std::size_t>(vehicle)];
            vehicle_arcs.reserve(
                static_cast<std::size_t>(location_count) *
                static_cast<std::size_t>(std::max(0, location_count - 1))
            );
            for (int tail = 0; tail < location_count; ++tail) {
                for (int head = 0; head < location_count; ++head) {
                    if (tail == head) {
                        continue;
                    }
                    VehicleAssignmentPhysicalArc arc;
                    arc.tail = tail;
                    arc.head = head;
                    arc.selected = model.addVar(
                        0.0,
                        1.0,
                        0.0,
                        GRB_BINARY,
                        "vi44_va_x_" + std::to_string(tail) + "_" +
                            std::to_string(head) + "_" + std::to_string(vehicle)
                    );
                    arc.flow = model.addVar(
                        0.0,
                        static_cast<double>(location_count),
                        0.0,
                        GRB_CONTINUOUS,
                        "vi44_va_f_" + std::to_string(tail) + "_" +
                            std::to_string(head) + "_" + std::to_string(vehicle)
                    );
                    vehicle_arcs.push_back(std::move(arc));
                }
            }
        }
    }

    model.update();
    GRBLinExpr feasibility_objective = 0.0;
    model.setObjective(feasibility_objective, GRB_MINIMIZE);

    for (std::size_t request_index = 0U; request_index < request_count; ++request_index) {
        GRBLinExpr pickup_cover = 0.0;
        GRBLinExpr delivery_cover = 0.0;
        for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
            pickup_cover +=
                pickup_assignment[request_index][static_cast<std::size_t>(vehicle)];
            delivery_cover +=
                delivery_assignment[request_index][static_cast<std::size_t>(vehicle)];
        }
        model.addConstr(
            pickup_cover == 1.0,
            "vi44_va_pickup_cover_" + std::to_string(request_index)
        );
        model.addConstr(
            delivery_cover == 1.0,
            "vi44_va_delivery_cover_" + std::to_string(request_index)
        );
    }

    for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
        GRBLinExpr pickup_count = 0.0;
        for (std::size_t request_index = 0U;
             request_index < request_count;
             ++request_index) {
            pickup_count +=
                pickup_assignment[request_index][static_cast<std::size_t>(vehicle)];
        }
        model.addConstr(
            pickup_count >= 1.0,
            "vi44_va_nonempty_" + std::to_string(vehicle)
        );
        if (vehicle + 1 < vehicle_count) {
            GRBLinExpr next_pickup_count = 0.0;
            for (std::size_t request_index = 0U;
                 request_index < request_count;
                 ++request_index) {
                next_pickup_count += pickup_assignment[request_index]
                    [static_cast<std::size_t>(vehicle + 1)];
            }
            model.addConstr(
                pickup_count >= next_pickup_count,
                "vi44_va_vehicle_size_order_" + std::to_string(vehicle)
            );
        }
    }
    for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
        for (std::size_t pickup_index = 0U;
             pickup_index < request_count;
             ++pickup_index) {
            GRBLinExpr matching_out = 0.0;
            for (const VehicleAssignmentMatchingArc& arc : matching_arcs) {
                if (arc.pickup_index == pickup_index) {
                    matching_out += arc.variables[static_cast<std::size_t>(vehicle)];
                }
            }
            model.addConstr(
                matching_out ==
                    pickup_assignment[pickup_index][static_cast<std::size_t>(vehicle)],
                "vi44_va_match_out_" + std::to_string(pickup_index) + "_" +
                    std::to_string(vehicle)
            );
        }
        for (std::size_t delivery_index = 0U;
             delivery_index < request_count;
             ++delivery_index) {
            GRBLinExpr matching_in = 0.0;
            for (const VehicleAssignmentMatchingArc& arc : matching_arcs) {
                if (arc.delivery_index == delivery_index) {
                    matching_in += arc.variables[static_cast<std::size_t>(vehicle)];
                }
            }
            model.addConstr(
                matching_in ==
                    delivery_assignment[delivery_index][static_cast<std::size_t>(vehicle)],
                "vi44_va_match_in_" + std::to_string(delivery_index) + "_" +
                    std::to_string(vehicle)
            );
        }
    }

    for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
        GRBLinExpr service_time = 0.0;
        for (std::size_t request_index = 0U;
             request_index < request_count;
             ++request_index) {
            service_time +=
                (data.time_pickup + data.time_empty) *
                pickup_assignment[request_index][static_cast<std::size_t>(vehicle)];
            service_time +=
                data.time_delivery *
                delivery_assignment[request_index][static_cast<std::size_t>(vehicle)];
        }

        if (options.add_container_bound) {
            GRBLinExpr container_workload = 0.0;
            for (std::size_t pickup_index = 0U;
                 pickup_index < request_count;
                 ++pickup_index) {
                const Request& pickup = data.requests[pickup_index];
                container_workload +=
                    data.time[static_cast<std::size_t>(pickup.from_id)]
                             [static_cast<std::size_t>(pickup.to_id)] *
                    pickup_assignment[pickup_index][static_cast<std::size_t>(vehicle)];
            }
            for (const VehicleAssignmentMatchingArc& arc : matching_arcs) {
                const Request& pickup = data.requests[arc.pickup_index];
                const Request& delivery = data.requests[arc.delivery_index];
                container_workload +=
                    data.time[static_cast<std::size_t>(pickup.to_id)]
                             [static_cast<std::size_t>(delivery.from_id)] *
                    arc.variables[static_cast<std::size_t>(vehicle)];
            }
            model.addConstr(
                service_time + 0.5 * container_workload <= data.time_limit,
                "vi44_va_container_time_" + std::to_string(vehicle)
            );
        }

        if (!options.add_tsp_bound) {
            continue;
        }

        for (int location = 0; location < location_count; ++location) {
            GRBLinExpr activity = 0.0;
            GRBLinExpr pickup_root_activity = 0.0;
            for (std::size_t request_index = 0U;
                 request_index < request_count;
                 ++request_index) {
                const Request& request = data.requests[request_index];
                if (request.from_id == location) {
                    activity += pickup_assignment[request_index]
                        [static_cast<std::size_t>(vehicle)];
                    activity += delivery_assignment[request_index]
                        [static_cast<std::size_t>(vehicle)];
                    pickup_root_activity += pickup_assignment[request_index]
                        [static_cast<std::size_t>(vehicle)];
                }
                if (request.to_id == location) {
                    activity += pickup_assignment[request_index]
                        [static_cast<std::size_t>(vehicle)];
                }
            }

            const GRBVar b =
                required_location[static_cast<std::size_t>(location)]
                                 [static_cast<std::size_t>(vehicle)];
            model.addConstr(
                b <= activity,
                "vi44_va_location_upper_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                path_start[static_cast<std::size_t>(location)]
                          [static_cast<std::size_t>(vehicle)] <= pickup_root_activity,
                "vi44_va_root_eligible_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                path_start[static_cast<std::size_t>(location)]
                          [static_cast<std::size_t>(vehicle)] <= b,
                "vi44_va_root_selected_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                path_end[static_cast<std::size_t>(location)]
                        [static_cast<std::size_t>(vehicle)] <= b,
                "vi44_va_end_selected_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );

            for (std::size_t request_index = 0U;
                 request_index < request_count;
                 ++request_index) {
                const Request& request = data.requests[request_index];
                if (request.from_id == location) {
                    model.addConstr(
                        b >= pickup_assignment[request_index]
                            [static_cast<std::size_t>(vehicle)],
                        "vi44_va_location_pickup_" + std::to_string(location) + "_" +
                            std::to_string(request_index) + "_" +
                            std::to_string(vehicle)
                    );
                    model.addConstr(
                        b >= delivery_assignment[request_index]
                            [static_cast<std::size_t>(vehicle)],
                        "vi44_va_location_delivery_" + std::to_string(location) + "_" +
                            std::to_string(request_index) + "_" +
                            std::to_string(vehicle)
                    );
                }
                if (request.to_id == location) {
                    model.addConstr(
                        b >= pickup_assignment[request_index]
                            [static_cast<std::size_t>(vehicle)],
                        "vi44_va_location_landfill_" + std::to_string(location) + "_" +
                            std::to_string(request_index) + "_" +
                            std::to_string(vehicle)
                    );
                }
            }
        }

        GRBLinExpr start_count = 0.0;
        GRBLinExpr end_count = 0.0;
        GRBLinExpr source_outflow = 0.0;
        GRBLinExpr selected_location_count = 0.0;
        for (int location = 0; location < location_count; ++location) {
            start_count += path_start[static_cast<std::size_t>(location)]
                                     [static_cast<std::size_t>(vehicle)];
            end_count += path_end[static_cast<std::size_t>(location)]
                                 [static_cast<std::size_t>(vehicle)];
            source_outflow += source_flow[static_cast<std::size_t>(location)]
                                         [static_cast<std::size_t>(vehicle)];
            selected_location_count +=
                required_location[static_cast<std::size_t>(location)]
                                 [static_cast<std::size_t>(vehicle)];
            model.addConstr(
                source_flow[static_cast<std::size_t>(location)]
                           [static_cast<std::size_t>(vehicle)] <=
                    static_cast<double>(location_count) *
                    path_start[static_cast<std::size_t>(location)]
                              [static_cast<std::size_t>(vehicle)],
                "vi44_va_source_flow_link_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
        }
        model.addConstr(start_count == 1.0, "vi44_va_one_start_" + std::to_string(vehicle));
        model.addConstr(end_count == 1.0, "vi44_va_one_end_" + std::to_string(vehicle));
        model.addConstr(
            source_outflow == selected_location_count,
            "vi44_va_source_balance_" + std::to_string(vehicle)
        );

        GRBLinExpr tsp_time = 0.0;
        auto& vehicle_arcs = physical_arcs[static_cast<std::size_t>(vehicle)];
        for (VehicleAssignmentPhysicalArc& arc : vehicle_arcs) {
            const GRBVar tail_selected =
                required_location[static_cast<std::size_t>(arc.tail)]
                                 [static_cast<std::size_t>(vehicle)];
            const GRBVar head_selected =
                required_location[static_cast<std::size_t>(arc.head)]
                                 [static_cast<std::size_t>(vehicle)];
            model.addConstr(
                arc.selected <= tail_selected,
                "vi44_va_arc_tail_" + std::to_string(arc.tail) + "_" +
                    std::to_string(arc.head) + "_" + std::to_string(vehicle)
            );
            model.addConstr(
                arc.selected <= head_selected,
                "vi44_va_arc_head_" + std::to_string(arc.tail) + "_" +
                    std::to_string(arc.head) + "_" + std::to_string(vehicle)
            );
            model.addConstr(
                arc.flow <= static_cast<double>(location_count) * arc.selected,
                "vi44_va_flow_arc_" + std::to_string(arc.tail) + "_" +
                    std::to_string(arc.head) + "_" + std::to_string(vehicle)
            );
            tsp_time +=
                data.time[static_cast<std::size_t>(arc.tail)]
                         [static_cast<std::size_t>(arc.head)] *
                arc.selected;
        }

        for (int location = 0; location < location_count; ++location) {
            GRBLinExpr incoming_arc = 0.0;
            GRBLinExpr outgoing_arc = 0.0;
            GRBLinExpr incoming_flow =
                source_flow[static_cast<std::size_t>(location)]
                           [static_cast<std::size_t>(vehicle)];
            GRBLinExpr outgoing_flow = 0.0;
            for (const VehicleAssignmentPhysicalArc& arc : vehicle_arcs) {
                if (arc.head == location) {
                    incoming_arc += arc.selected;
                    incoming_flow += arc.flow;
                }
                if (arc.tail == location) {
                    outgoing_arc += arc.selected;
                    outgoing_flow += arc.flow;
                }
            }
            const GRBVar b =
                required_location[static_cast<std::size_t>(location)]
                                 [static_cast<std::size_t>(vehicle)];
            model.addConstr(
                path_start[static_cast<std::size_t>(location)]
                          [static_cast<std::size_t>(vehicle)] +
                    incoming_arc == b,
                "vi44_va_path_in_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                path_end[static_cast<std::size_t>(location)]
                        [static_cast<std::size_t>(vehicle)] +
                    outgoing_arc == b,
                "vi44_va_path_out_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
            model.addConstr(
                incoming_flow - outgoing_flow == b,
                "vi44_va_flow_balance_" + std::to_string(location) + "_" +
                    std::to_string(vehicle)
            );
        }
        model.addConstr(
            service_time + tsp_time <= data.time_limit,
            "vi44_va_tsp_time_" + std::to_string(vehicle)
        );
    }

    const auto solve_start = std::chrono::steady_clock::now();
    model.optimize();
    const auto solve_end = std::chrono::steady_clock::now();

    VehicleAssignmentFeasibilityResult result;
    result.status = model.get(GRB_IntAttr_Status);
    result.hit_time_limit = result.status == GRB_TIME_LIMIT;
    result.runtime_seconds =
        std::chrono::duration<double>(solve_end - solve_start).count();
    result.feasible = model.get(GRB_IntAttr_SolCount) > 0;
    result.proven_infeasible = result.status == GRB_INFEASIBLE;

    if (!result.feasible && !result.proven_infeasible &&
        result.status != GRB_TIME_LIMIT) {
        throw std::runtime_error(
            "The VI-44 vehicle-assignment problem stopped with unsupported status " +
            std::to_string(result.status) + "."
        );
    }

    if (!result.feasible) {
        return result;
    }

    result.vehicles.reserve(static_cast<std::size_t>(vehicle_count));
    for (int vehicle = 0; vehicle < vehicle_count; ++vehicle) {
        VI44VehicleAssignmentVehicleResult vehicle_result;
        vehicle_result.vehicle_index = vehicle;
        for (std::size_t request_index = 0U;
             request_index < request_count;
             ++request_index) {
            if (pickup_assignment[request_index][static_cast<std::size_t>(vehicle)]
                    .get(GRB_DoubleAttr_X) > 0.5) {
                ++vehicle_result.pickup_count;
            }
            if (delivery_assignment[request_index][static_cast<std::size_t>(vehicle)]
                    .get(GRB_DoubleAttr_X) > 0.5) {
                ++vehicle_result.delivery_count;
            }
        }
        vehicle_result.service_time =
            (data.time_pickup + data.time_empty) *
                static_cast<double>(vehicle_result.pickup_count) +
            data.time_delivery * static_cast<double>(vehicle_result.delivery_count);

        double active_travel_lb = 0.0;
        if (options.add_container_bound) {
            double workload = 0.0;
            for (std::size_t pickup_index = 0U;
                 pickup_index < request_count;
                 ++pickup_index) {
                if (pickup_assignment[pickup_index][static_cast<std::size_t>(vehicle)]
                        .get(GRB_DoubleAttr_X) <= 0.5) {
                    continue;
                }
                const Request& pickup = data.requests[pickup_index];
                workload +=
                    data.time[static_cast<std::size_t>(pickup.from_id)]
                             [static_cast<std::size_t>(pickup.to_id)];
            }
            for (const VehicleAssignmentMatchingArc& arc : matching_arcs) {
                if (arc.variables[static_cast<std::size_t>(vehicle)]
                        .get(GRB_DoubleAttr_X) <= 0.5) {
                    continue;
                }
                const Request& pickup = data.requests[arc.pickup_index];
                const Request& delivery = data.requests[arc.delivery_index];
                workload +=
                    data.time[static_cast<std::size_t>(pickup.to_id)]
                             [static_cast<std::size_t>(delivery.from_id)];
            }
            vehicle_result.container_lb = 0.5 * workload;
            active_travel_lb = vehicle_result.container_lb;
        }
        if (options.add_tsp_bound) {
            double tsp_lb = 0.0;
            for (const VehicleAssignmentPhysicalArc& arc :
                 physical_arcs[static_cast<std::size_t>(vehicle)]) {
                if (arc.selected.get(GRB_DoubleAttr_X) > 0.5) {
                    tsp_lb +=
                        data.time[static_cast<std::size_t>(arc.tail)]
                                 [static_cast<std::size_t>(arc.head)];
                }
            }
            vehicle_result.tsp_lb = tsp_lb;
            active_travel_lb = std::max(active_travel_lb, tsp_lb);
        }
        vehicle_result.active_duration_lb =
            vehicle_result.service_time + active_travel_lb;
        result.vehicles.push_back(vehicle_result);
    }
    return result;
}

VI44KMinResult compute_vehicle_assignment_k_min(
    const SPDPData& data,
    const VI44VehicleAssignmentOptions& options
) {
    if (!options.add_tsp_bound && !options.add_container_bound) {
        throw std::runtime_error(
            "VI-44 vehicle assignment requires the TSP bound, the container bound, "
            "or both."
        );
    }
    if (!std::isfinite(options.time_limit) || options.time_limit < 0.0) {
        throw std::runtime_error(
            "The VI-44 vehicle-assignment time limit must be finite and nonnegative."
        );
    }

    VI44KMinResult result;
    result.vehicle_assignment_enabled = true;
    if (data.requests.empty()) {
        result.vehicle_assignment_has_certified_bound = true;
        result.vehicle_assignment_optimal = true;
        return result;
    }

    GRBEnv environment(true);
    environment.set(GRB_IntParam_OutputFlag, 0);
    environment.start();

    const auto overall_start = std::chrono::steady_clock::now();
    int certified_k_min = 1;
    for (int vehicle_count = 1;
         vehicle_count <= static_cast<int>(data.requests.size());
         ++vehicle_count) {
        double remaining_time = 0.0;
        if (options.time_limit > 0.0) {
            const double elapsed =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - overall_start
                ).count();
            remaining_time = std::max(0.0, options.time_limit - elapsed);
            if (remaining_time <= kTolerance) {
                result.vehicle_assignment_hit_time_limit = true;
                result.vehicle_assignment_k_min = certified_k_min;
                result.vehicle_assignment_has_certified_bound = true;
                result.vehicle_assignment_runtime_seconds = elapsed;
                return result;
            }
        }

        const VehicleAssignmentFeasibilityResult feasibility =
            solve_vehicle_assignment_feasibility(
                data,
                vehicle_count,
                options,
                environment,
                remaining_time
            );
        result.vehicle_assignment_status = feasibility.status;
        result.vehicle_assignment_tested_vehicle_count = vehicle_count;

        if (feasibility.feasible) {
            result.vehicle_assignment_k_min = vehicle_count;
            result.vehicle_assignment_has_certified_bound = true;
            result.vehicle_assignment_has_feasible_solution = true;
            result.vehicle_assignment_optimal = true;
            result.vehicle_assignment_vehicles = feasibility.vehicles;
            break;
        }
        if (feasibility.proven_infeasible) {
            certified_k_min = vehicle_count + 1;
            continue;
        }

        result.vehicle_assignment_hit_time_limit = feasibility.hit_time_limit;
        result.vehicle_assignment_k_min = certified_k_min;
        result.vehicle_assignment_has_certified_bound = true;
        break;
    }

    result.vehicle_assignment_runtime_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - overall_start
        ).count();
    if (!result.vehicle_assignment_has_certified_bound) {
        result.vehicle_assignment_k_min = certified_k_min;
        result.vehicle_assignment_has_certified_bound = true;
    }
    return result;
}

}  // namespace

VI44KMinResult compute_vi44_k_min(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const VI44KMinOptions& options
) {
    if (options.precomputed_result != nullptr) {
        return *options.precomputed_result;
    }
    if (!options.use_cor && !options.use_subproblem &&
        !options.use_vehicle_assignment) {
        throw std::runtime_error(
            "VI-44 is enabled, but every k_min method is disabled."
        );
    }

    VI44KMinResult result;
    if (options.use_cor) {
        result.cor_enabled = true;
        result.cor_k_min = cor_route_lower_bound_kmin(data);
        result.selected_k_min = std::max(result.selected_k_min, result.cor_k_min);
    }

    if (options.use_subproblem) {
        if (!std::isfinite(options.subproblem.time_limit) ||
            options.subproblem.time_limit < 0.0) {
            throw std::runtime_error(
                "The auxiliary VI-44 subproblem time limit must be finite and "
                "nonnegative."
            );
        }
        result.subproblem_enabled = true;
        const AuxiliaryDurationSubproblemResult auxiliary =
            solve_auxiliary_duration_subproblem(
                data,
                graph,
                options.subproblem.time_limit,
                options.subproblem.type == VI44SubproblemType::IP,
                options.subproblem.add_time_constraints,
                options.subproblem.rounded_bound_stop
            );
        result.subproblem_status = auxiliary.status;
        result.subproblem_hit_time_limit = auxiliary.hit_time_limit;
        result.subproblem_has_certified_bound = auxiliary.has_certified_bound;
        result.subproblem_objective_value = auxiliary.objective_value;
        result.subproblem_objective_bound = auxiliary.objective_bound;
        result.subproblem_runtime_seconds = auxiliary.runtime_seconds;
        result.subproblem_stopped_by_rounded_bound =
            auxiliary.stopped_by_rounded_bound;
        result.subproblem_rounded_bound_certified =
            auxiliary.rounded_bound_certified;
        result.subproblem_certified_rounded_k =
            auxiliary.certified_rounded_k;
        result.subproblem_callback_objective_ub =
            auxiliary.callback_objective_ub;
        result.subproblem_callback_safe_objective_lb =
            auxiliary.callback_safe_objective_lb;

        if (auxiliary.has_certified_bound) {
            const double scale = std::max(
                {1.0, data.time_limit, std::abs(auxiliary.objective_bound)}
            );
            result.subproblem_numerical_tolerance =
                10.0 * kAuxiliaryLPSolverTolerance * scale;
            result.subproblem_safe_lower_bound = std::max(
                0.0,
                auxiliary.objective_bound -
                    result.subproblem_numerical_tolerance
            );
            result.subproblem_k_min = std::max(
                0,
                static_cast<int>(
                    std::ceil(
                        result.subproblem_safe_lower_bound / data.time_limit
                    )
                )
            );
            result.selected_k_min =
                std::max(result.selected_k_min, result.subproblem_k_min);
        }
    }

    if (options.use_vehicle_assignment) {
        VI44KMinResult vehicle_assignment_result =
            compute_vehicle_assignment_k_min(data, options.vehicle_assignment);
        result.vehicle_assignment_enabled = true;
        result.vehicle_assignment_k_min =
            vehicle_assignment_result.vehicle_assignment_k_min;
        result.vehicle_assignment_status =
            vehicle_assignment_result.vehicle_assignment_status;
        result.vehicle_assignment_hit_time_limit =
            vehicle_assignment_result.vehicle_assignment_hit_time_limit;
        result.vehicle_assignment_has_certified_bound =
            vehicle_assignment_result.vehicle_assignment_has_certified_bound;
        result.vehicle_assignment_has_feasible_solution =
            vehicle_assignment_result.vehicle_assignment_has_feasible_solution;
        result.vehicle_assignment_optimal =
            vehicle_assignment_result.vehicle_assignment_optimal;
        result.vehicle_assignment_tested_vehicle_count =
            vehicle_assignment_result.vehicle_assignment_tested_vehicle_count;
        result.vehicle_assignment_runtime_seconds =
            vehicle_assignment_result.vehicle_assignment_runtime_seconds;
        result.vehicle_assignment_vehicles =
            std::move(vehicle_assignment_result.vehicle_assignment_vehicles);
        if (result.vehicle_assignment_has_certified_bound) {
            result.selected_k_min =
                std::max(result.selected_k_min, result.vehicle_assignment_k_min);
        }
    }

    return result;
}

std::optional<int> certified_rounded_duration_bound(
    double exact_incumbent_upper_bound,
    double solver_objective_lower_bound,
    double route_time_limit
) {
    if (!std::isfinite(exact_incumbent_upper_bound) ||
        !std::isfinite(solver_objective_lower_bound) ||
        !std::isfinite(route_time_limit) || route_time_limit <= 0.0 ||
        exact_incumbent_upper_bound < 0.0) {
        return std::nullopt;
    }
    const double safe_lower = safe_duration_lower_bound(
        solver_objective_lower_bound, route_time_limit);
    const int lower_k = static_cast<int>(std::ceil(safe_lower / route_time_limit));
    const int upper_k = static_cast<int>(
        std::ceil(exact_incumbent_upper_bound / route_time_limit));
    if (lower_k != upper_k) {
        return std::nullopt;
    }
    return lower_k;
}

}  // namespace spdp
