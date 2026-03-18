#include "PstepFormulation.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ios>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace spdp {
namespace {

constexpr double kTolerance = 1e-9;

// 부동소수점 비교에서 사용하는 허용 오차 기반 동등 비교.
bool double_equal(double lhs, double rhs) {
    return std::fabs(lhs - rhs) <= kTolerance;
}

// lhs <= rhs 를 허용 오차와 함께 판정한다.
bool double_less_or_equal(double lhs, double rhs) {
    return lhs <= rhs + kTolerance;
}

// lhs >= rhs 를 허용 오차와 함께 판정한다.
bool double_greater_or_equal(double lhs, double rhs) {
    return lhs + kTolerance >= rhs;
}

// 상태 슬롯 순서를 정규화해 같은 의미의 상태가 동일하게 비교되도록 맞춘다.
State canonicalize_state(State state) {
    if (state[1] < state[0]) {
        std::swap(state[0], state[1]);
    }
    return state;
}

// 정규화된 상태 기준 사전식 비교.
bool state_less(const State& lhs, const State& rhs) {
    const State lhs_canonical = canonicalize_state(lhs);
    const State rhs_canonical = canonicalize_state(rhs);

    if (lhs_canonical[0] < rhs_canonical[0]) {
        return true;
    }
    if (rhs_canonical[0] < lhs_canonical[0]) {
        return false;
    }
    return lhs_canonical[1] < rhs_canonical[1];
}

// 슬롯 순서와 무관하게 두 상태가 같은지 판정한다.
bool same_state(const State& lhs, const State& rhs) {
    return canonicalize_state(lhs) == canonicalize_state(rhs);
}

// 로그 출력용 실수 포맷팅.
std::string format_double(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

// 조건이 깨지면 일관된 방식으로 예외를 던진다.
void require_condition(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

// raw p-step 길이 조건을 만족하는지 확인한다.
bool is_valid_raw_length(
    NodeId start_node_id,
    NodeId last_node_id,
    int q,
    int p,
    NodeId end_node_id
) {
    if (start_node_id == 0) {
        return q >= 1 && q <= p;
    }
    return q == p;
}

// 현재 경로를 raw p-step으로 기록해야 하는지 결정한다.
bool should_record_raw_path(
    NodeId start_node_id,
    NodeId last_node_id,
    int q,
    int p,
    NodeId end_node_id
) {
    return is_valid_raw_length(start_node_id, last_node_id, q, p, end_node_id);
}

// edge id 시퀀스로부터 RawPStepPath를 구성한다.
RawPStepPath make_raw_path(
    int path_id,
    const std::vector<int>& edge_ids,
    const std::vector<EdgeRecord>& all_edges
) {
    require_condition(!edge_ids.empty(), "Raw p-step path must contain at least one edge.");

    RawPStepPath path;
    path.id = path_id;
    path.edge_ids = edge_ids;
    path.q = static_cast<int>(edge_ids.size());

    const EdgeRecord& first_edge = all_edges[static_cast<std::size_t>(edge_ids.front())];
    const EdgeRecord& last_edge = all_edges[static_cast<std::size_t>(edge_ids.back())];

    path.start_node_id = first_edge.u;
    path.last_node_id = last_edge.v;
    path.start_state = canonicalize_state(first_edge.data.start_state);
    path.last_state = canonicalize_state(last_edge.data.end_state);
    path.node_sequence.reserve(edge_ids.size() + 1U);
    path.state_sequence.reserve(edge_ids.size() + 1U);
    path.node_sequence.push_back(first_edge.u);
    path.state_sequence.push_back(path.start_state);

    for (int edge_id : edge_ids) {
        const EdgeRecord& edge = all_edges[static_cast<std::size_t>(edge_id)];
        path.total_time += edge.data.time;
        path.total_cost += edge.data.cost;
        path.node_sequence.push_back(edge.v);
        path.state_sequence.push_back(canonicalize_state(edge.data.end_state));
    }

    return path;
}

// 하나의 raw path에서 허용되는 극단 tau 값들만 추린다.
std::vector<double> extreme_tau_values(
    const RawPStepPath& path,
    double time_limit,
    NodeId end_node_id
) {
    std::vector<double> tau_values;
    const double latest_start = time_limit - path.total_time;

    if (path.start_node_id == 0) {
        tau_values.push_back(0.0);
    } else if (path.last_node_id == end_node_id) {
        if (double_greater_or_equal(latest_start, 0.0)) {
            tau_values.push_back(std::max(0.0, latest_start));
        }
    } else {
        tau_values.push_back(0.0);
        if (double_greater_or_equal(latest_start, 0.0)) {
            tau_values.push_back(std::max(0.0, latest_start));
        }
    }

    std::sort(tau_values.begin(), tau_values.end());
    tau_values.erase(
        std::unique(
            tau_values.begin(),
            tau_values.end(),
            [](double lhs, double rhs) {
                return double_equal(lhs, rhs);
            }
        ),
        tau_values.end()
    );
    return tau_values;
}

// edge id 목록을 디버깅용 문자열로 바꾼다.
std::string format_edge_sequence(const std::vector<int>& edge_ids) {
    std::ostringstream out;
    out << "[";
    for (std::size_t idx = 0; idx < edge_ids.size(); ++idx) {
        if (idx > 0U) {
            out << ", ";
        }
        out << edge_ids[idx];
    }
    out << "]";
    return out.str();
}

// 활성 theta edge를 사람이 읽기 쉬운 디버그 문자열로 만든다.
std::string format_debug_edge(
    int edge_id,
    const EdgeRecord& edge
) {
    std::ostringstream out;
    out << "theta_" << edge_id
        << ": " << edge.u << "->" << edge.v
        << " pi=" << format_edge_sequence(edge.data.sequence_pi)
        << " start=" << state_to_str(edge.data.start_state)
        << " end=" << state_to_str(edge.data.end_state);
    return out.str();
}

// 양의 값을 갖는 x 변수용 디버그 문자열을 만든다.
std::string format_debug_pstep(
    const CompactPStep& pstep,
    double value
) {
    std::ostringstream out;
    out << "x_" << pstep.id
        << "=" << format_double(value)
        << " tau=" << format_double(pstep.tau)
        << " edges=" << format_edge_sequence(pstep.edge_ids)
        << " start=" << state_to_str(pstep.start_state)
        << " end=" << state_to_str(pstep.last_state);
    return out.str();
}

// node 시퀀스를 출력용 문자열로 변환한다.
std::string format_node_sequence(const std::vector<NodeId>& node_sequence) {
    std::ostringstream out;
    out << "[";
    for (std::size_t idx = 0; idx < node_sequence.size(); ++idx) {
        if (idx > 0U) {
            out << ", ";
        }
        out << node_sequence[idx];
    }
    out << "]";
    return out.str();
}

// 상태 전이 시퀀스를 출력용 문자열로 변환한다.
std::string format_state_sequence(const std::vector<State>& state_sequence) {
    std::ostringstream out;
    out << "[";
    for (std::size_t idx = 0; idx < state_sequence.size(); ++idx) {
        if (idx > 0U) {
            out << " -> ";
        }
        out << state_to_str(state_sequence[idx]);
    }
    out << "]";
    return out.str();
}

// 방문 제약식에 들어갈 계수를 p-step에서 계산한다.
std::map<NodeId, int> calculated_visit_coefficients(
    const CompactPStep& pstep,
    const MultiDiGraph& graph
) {
    std::map<NodeId, int> coefficients;
    for (std::size_t pos = 0; pos < pstep.node_sequence.size(); ++pos) {
        const NodeId node_id = pstep.node_sequence[pos];
        if (!graph.is_physical_service_node(node_id)) {
            continue;
        }

        coefficients[node_id] = (pos == 0U || pos + 1U == pstep.node_sequence.size()) ? 1 : 2;
    }
    return coefficients;
}

// 상태 보존 제약식에 들어갈 계수를 p-step에서 계산한다.
std::map<NodeStateKey, int> calculated_state_coefficients(
    const CompactPStep& pstep,
    const MultiDiGraph& graph
) {
    std::map<NodeStateKey, int> coefficients;

    if (graph.is_physical_service_node(pstep.start_node_id)) {
        coefficients[{pstep.start_node_id, canonicalize_state(pstep.start_state)}] = 1;
    }
    if (graph.is_physical_service_node(pstep.last_node_id)) {
        coefficients[{pstep.last_node_id, canonicalize_state(pstep.last_state)}] = -1;
    }

    return coefficients;
}

// 시간 관련 제약식에 들어갈 계수를 p-step에서 계산한다.
std::map<NodeStateKey, double> calculated_time_coefficients(
    const CompactPStep& pstep,
    const MultiDiGraph& graph
) {
    std::map<NodeStateKey, double> coefficients;

    if (graph.is_physical_service_node(pstep.start_node_id)) {
        coefficients[{pstep.start_node_id, canonicalize_state(pstep.start_state)}] = pstep.tau;
    }
    if (graph.is_physical_service_node(pstep.last_node_id)) {
        coefficients[{pstep.last_node_id, canonicalize_state(pstep.last_state)}] =
            -(pstep.tau + pstep.total_time);
    }

    return coefficients;
}

// 해 복원 과정에서 차량에 실린 요청 상태를 추적하기 위한 내부 구조체.
struct OnboardRequest {
    int request_id = 0;
    int container_type = -1;
    int treatment_location = 0;
    bool is_full = true;
};

// 해 출력 시 정수값처럼 보이는 실수는 깔끔하게 정리해 출력한다.
std::string format_solution_value(double value) {
    const double rounded = std::round(value);
    if (double_equal(value, rounded)) {
        return std::to_string(static_cast<long long>(rounded));
    }
    return format_double(value);
}

// theta 변수값이 1인 edge id만 모은다.
std::vector<int> collect_active_edge_ids(const CompactMasterProblem& problem) {
    std::vector<int> active_edge_ids;
    for (std::size_t edge_id = 0; edge_id < problem.theta_vars.size(); ++edge_id) {
        if (problem.theta_vars[edge_id].get(GRB_DoubleAttr_X) > 0.5) {
            active_edge_ids.push_back(static_cast<int>(edge_id));
        }
    }
    return active_edge_ids;
}

// 현재 적재 상태를 onboard 요청 목록으로부터 State 형태로 변환한다.
// 해 복원 중에는 "차량에 어떤 요청들이 실려 있는가"를 직접 추적하고,
// 그 결과가 multigraph edge에 저장된 start/end state와 일치하는지 계속 확인한다.
// 이 함수는 그 연결 고리 역할을 한다.
State make_state_from_onboard(
    const std::vector<OnboardRequest>& onboard_requests
) {
    require_condition(
        onboard_requests.size() <= 2U,
        "Recovered onboard request count exceeds vehicle capacity."
    );

    std::vector<StateToken> tokens;
    tokens.reserve(2);

    for (const OnboardRequest& request : onboard_requests) {
        if (request.is_full) {
            tokens.push_back(StateToken{'F', request.container_type, request.treatment_location});
        } else {
            tokens.push_back(StateToken{'E', request.container_type, -1});
        }
    }

    while (tokens.size() < 2U) {
        tokens.push_back(StateToken{'N', -1, -1});
    }

    State state{tokens[0], tokens[1]};
    return canonicalize_state(state);
}

// 하나의 edge를 현재 onboard 상태에 적용할 수 있는지 확인하고, 가능하면 다음 상태를 만든다.
// 즉 "현재 차량 적재 상태에서 이 edge를 실제로 탈 수 있는가"를 시뮬레이션하는 함수다.
// 성공하면 next_onboard는 edge 적용 후의 적재 상태가 된다.
bool try_apply_edge_to_onboard(
    const EdgeRecord& edge,
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<OnboardRequest>& current_onboard,
    std::vector<OnboardRequest>& next_onboard
) {
    const State current_state = make_state_from_onboard(current_onboard);
    if (!same_state(current_state, edge.data.start_state)) {
        return false;
    }

    next_onboard = current_onboard;

    // edge 내부의 treatment sequence를 따라 full container를 empty로 바꾼다.
    for (int treatment_location : edge.data.sequence_pi) {
        for (auto it = next_onboard.rbegin(); it != next_onboard.rend(); ++it) {
            if (it->is_full && it->treatment_location == treatment_location) {
                it->is_full = false;
            }
        }
    }

    const NodeSpec& arrival_node = graph.node(edge.v);
    if (arrival_node.kind == NodeSpec::Kind::Pickup) {
        // Pickup node에 도착하면 새로운 full container 요청이 onboard에 추가된다.
        next_onboard.push_back(OnboardRequest{
            arrival_node.request_idx.value(),
            arrival_node.container_type.value(),
            arrival_node.landfill_location.value(),
            true,
        });
    } else if (arrival_node.kind == NodeSpec::Kind::Delivery) {
        // Delivery node에 도착하면 대응되는 empty container 하나를 내려야 한다.
        const int container_type = arrival_node.container_type.value();
        auto found = next_onboard.end();
        for (auto it = next_onboard.begin(); it != next_onboard.end(); ++it) {
            if (!it->is_full && it->container_type == container_type) {
                found = it;
            }
        }
        if (found == next_onboard.end()) {
            return false;
        }
        next_onboard.erase(found);
    }

    // 시뮬레이션 결과가 edge가 요구하는 end_state와 일치해야만 이 edge를 탈 수 있다.
    const State next_state = make_state_from_onboard(next_onboard);
    if (!same_state(next_state, edge.data.end_state)) {
        return false;
    }

    return true;
}

bool decompose_selected_edges(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<EdgeRecord>& all_edges,
    const std::map<NodeId, std::vector<int>>& active_outgoing_edges,
    std::vector<bool>& is_edge_consumed,
    std::vector<std::vector<int>>& route_edge_sets,
    NodeId end_node_id
);

bool extend_current_route(
    NodeId current_node_id,
    const std::vector<OnboardRequest>& current_onboard,
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<EdgeRecord>& all_edges,
    const std::map<NodeId, std::vector<int>>& active_outgoing_edges,
    std::vector<bool>& is_edge_consumed,
    std::vector<int>& current_route_edge_ids,
    std::vector<std::vector<int>>& route_edge_sets,
    NodeId end_node_id
) {
    // 현재까지 선택한 active theta edges를 따라 한 개의 완전한 route를 만들려고 시도한다.
    // route가 하나 완성되면 남은 active edge들도 같은 방식으로 분해 가능한지 decompose_selected_edges를 호출하여 재귀적으로 확인한다.
    // 반환값의 의미:
    // - true  = 현재까지의 선택을 포함해, 남은 active edge 전체를 끝까지 route들로 분해할 수 있음
    // - false = 현재 선택으로는 전체 분해가 불가능하므로, 이전 단계로 돌아가 다른 edge를 시도해야 함

    // 종료 노드에 도달했고 차량이 비어 있으면 하나의 route 분해가 완성된다.
    if (current_node_id == end_node_id) {
        if (!current_onboard.empty()) {
            return false;
        }

        // current_route_edge_ids 하나가 완성되었으므로 route 집합에 잠시 추가한다.
        // 그 다음 남은 미사용 edge들도 전부 분해 가능한지 decompose_selected_edges()로 다시 확인한다.
        route_edge_sets.push_back(current_route_edge_ids);
        if (decompose_selected_edges(
                data,
                graph,
                all_edges,
                active_outgoing_edges,
                is_edge_consumed,
                route_edge_sets,
                end_node_id
            )) {
            return true;
        }
        route_edge_sets.pop_back();
        return false;
    }

    const auto found = active_outgoing_edges.find(current_node_id);
    if (found == active_outgoing_edges.end()) {
        return false;
    }

    std::vector<int> candidate_edge_ids;
    const State current_state = make_state_from_onboard(current_onboard);
    // 현재 노드에서 나가는 active edge 중에서도,
    // onboard 상태가 edge.start_state와 맞는 것만 실제 후보가 된다.
    for (int edge_id : found->second) {
        if (is_edge_consumed[static_cast<std::size_t>(edge_id)]) {
            continue;
        }
        if (same_state(all_edges[static_cast<std::size_t>(edge_id)].data.start_state, current_state)) {
            candidate_edge_ids.push_back(edge_id);
        }
    }

    std::sort(candidate_edge_ids.begin(), candidate_edge_ids.end());

    // 후보 edge를 하나 고른 뒤 실제 onboard 상태 전이를 적용해 보고,
    // 끝까지 모순 없이 이어지면 그 선택을 채택한다.
    // 중간에 막히면 백트래킹으로 이전 상태로 되돌아간다.
    for (int edge_id : candidate_edge_ids) {
        std::vector<OnboardRequest> next_onboard;
        if (!try_apply_edge_to_onboard(
                all_edges[static_cast<std::size_t>(edge_id)],
                data,
                graph,
                current_onboard,
                next_onboard
            )) {
            continue;
        }

        is_edge_consumed[static_cast<std::size_t>(edge_id)] = true;
        current_route_edge_ids.push_back(edge_id);

        const NodeId next_node_id = all_edges[static_cast<std::size_t>(edge_id)].v;
        // 현재 route를 한 edge 더 연장한 상태로 재귀 호출한다.
        // 여기서 true가 돌아오면 "이 edge를 고른 선택이 전체 분해 성공으로 이어진다"는 뜻이다.
        if (extend_current_route(
                next_node_id,
                next_onboard,
                data,
                graph,
                all_edges,
                active_outgoing_edges,
                is_edge_consumed,
                current_route_edge_ids,
                route_edge_sets,
                end_node_id
            )) {
            return true;
        }

        // 이 edge를 고른 경우 전체 분해가 실패했으므로, 상태를 원래대로 되돌리고 다음 후보를 본다.
        current_route_edge_ids.pop_back();
        is_edge_consumed[static_cast<std::size_t>(edge_id)] = false;
    }

    // 어떤 후보 edge를 선택해도 전체 분해가 완성되지 않으면 실패를 반환한다.
    return false;
}

bool decompose_selected_edges(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<EdgeRecord>& all_edges,
    const std::map<NodeId, std::vector<int>>& active_outgoing_edges,
    std::vector<bool>& is_edge_consumed,
    std::vector<std::vector<int>>& route_edge_sets,
    NodeId end_node_id
) {
    // active theta edge 전체를 여러 개의 feasible route로 분해할 수 있는지 검사한다.
    // 해 복원에서는 이 함수가 핵심이며, edge support가 실제 route 집합에 대응하는지 확인한다.
    // 반환값의 의미:
    // - true  = 현재 남아 있는 미사용 edge들을 전부 route들로 분해하는 데 성공
    // - false = 어떤 식으로 route를 나눠도 전부 소진하는 데 실패
    
    // 더 이상 미사용 edge가 없으면 분해 성공이다.
    bool has_unused_edge = false;
    for (bool consumed : is_edge_consumed) {
        if (!consumed) {
            has_unused_edge = true;
            break;
        }
    }
    if (!has_unused_edge) {
        return true;
    }

    // 아직 미사용 edge가 남아 있으면 depot에서 새 route 하나를 시작해 본다.
    // 이 route가 완성되면 extend_current_route() 안에서 다시 decompose_selected_edges()가 호출되어
    // "남은 edge들"의 분해 가능 여부를 재귀적으로 확인한다.
    std::vector<int> current_route_edge_ids;
    return extend_current_route(
        0,
        {},
        data,
        graph,
        all_edges,
        active_outgoing_edges,
        is_edge_consumed,
        current_route_edge_ids,
        route_edge_sets,
        end_node_id
    );
}

// 하나의 edge 시퀀스를 사람이 읽을 수 있는 route 정보로 복원한다.
// edge 수준의 경로를 Path / Loc / Patt / Load 형식의 최종 출력용 route action들로 바꾼다.
RecoveredRouteSolution build_route_solution(
    int route_id,
    const std::vector<int>& route_edge_ids,
    const SPDPData& data,
    const MultiDiGraph& graph
) {
    const std::vector<EdgeRecord>& all_edges = graph.edges();
    const NodeId end_node_id = graph.end_node_id();
    require_condition(!route_edge_ids.empty(), "Recovered route must contain at least one edge.");

    RecoveredRouteSolution route;
    route.route_id = route_id;
    route.edge_ids = route_edge_ids;

    std::vector<OnboardRequest> onboard_requests;
    NodeId current_node_id = 0;
    int load = 0;

    const EdgeRecord& first_edge = all_edges[static_cast<std::size_t>(route_edge_ids.front())];
    require_condition(first_edge.u == 0, "Recovered route does not start at the start depot.");

    for (int edge_id : route_edge_ids) {
        require_condition(edge_id >= 0, "Recovered route edge id must be nonnegative.");
        require_condition(
            static_cast<std::size_t>(edge_id) < all_edges.size(),
            "Recovered route edge id is out of range."
        );

        const EdgeRecord& edge = all_edges[static_cast<std::size_t>(edge_id)];
        require_condition(edge.u == current_node_id, "Recovered route violates node continuity.");

        route.total_time += edge.data.time;
        route.total_cost += edge.data.cost;

        // edge 내부 sequence_pi에 따라 full container가 treatment를 거쳐 empty가 된다.
        // 출력에서는 pattern=1로 기록한다.
        for (int treatment_location : edge.data.sequence_pi) {
            bool emptied_any_container = false;
            for (auto it = onboard_requests.rbegin(); it != onboard_requests.rend(); ++it) {
                if (!it->is_full || it->treatment_location != treatment_location) {
                    continue;
                }

                emptied_any_container = true;
                route.actions.push_back(RouteSolutionAction{
                    it->request_id,
                    treatment_location,
                    1,
                    load,
                });
                it->is_full = false;
            }

            require_condition(
                emptied_any_container,
                "Recovered route visits a treatment location without a matching full container onboard."
            );
        }

        const NodeSpec& arrival_node = graph.node(edge.v);
        if (arrival_node.kind == NodeSpec::Kind::Pickup) {
            // Pickup 도착 시 새 요청이 차량에 적재된다.
            // 출력에서는 pattern=0으로 기록한다.
            require_condition(
                onboard_requests.size() < 2U,
                "Recovered route exceeds vehicle capacity at pickup."
            );
            const int request_id = arrival_node.request_idx.value();
            onboard_requests.push_back(OnboardRequest{
                request_id,
                arrival_node.container_type.value(),
                arrival_node.landfill_location.value(),
                true,
            });
            ++load;
            route.actions.push_back(RouteSolutionAction{
                request_id,
                arrival_node.location,
                0,
                load,
            });
        } else if (arrival_node.kind == NodeSpec::Kind::Delivery) {
            // Delivery 도착 시 대응되는 empty container 요청을 차량에서 내린다.
            // 출력에서는 pattern=2로 기록한다.
            const int request_id = arrival_node.request_idx.value();
            const int container_type = arrival_node.container_type.value();
            auto found = onboard_requests.end();
            for (auto it = onboard_requests.begin(); it != onboard_requests.end(); ++it) {
                if (!it->is_full && it->container_type == container_type) {
                    found = it;
                }
            }

            require_condition(
                found != onboard_requests.end(),
                "Recovered route tries to deliver without a matching empty container onboard."
            );

            onboard_requests.erase(found);
            --load;
            route.actions.push_back(RouteSolutionAction{
                request_id,
                arrival_node.location,
                2,
                load,
            });
        }

        require_condition(
            load >= 0 && load <= 2,
            "Recovered route violates the vehicle load bounds."
        );
        require_condition(
            load == static_cast<int>(onboard_requests.size()),
            "Recovered route load is inconsistent with the onboard container count."
        );
        current_node_id = edge.v;
    }

    require_condition(
        current_node_id == end_node_id,
        "Recovered route does not end at the end depot."
    );
    require_condition(
        double_less_or_equal(route.total_time, data.time_limit),
        "Recovered route exceeds the instance time limit."
    );
    require_condition(
        onboard_requests.empty(),
        "Recovered route ends with containers still onboard."
    );
    require_condition(load == 0, "Recovered route ended with a nonzero onboard load.");
    require_condition(
        std::all_of(
            onboard_requests.begin(),
            onboard_requests.end(),
            [](const OnboardRequest& request) {
                return !request.is_full;
            }
        ),
        "Recovered route still contains full requests at the end."
    );

    return route;
}

}  // namespace

bool NodeStateKey::operator<(const NodeStateKey& other) const {
    if (node_id != other.node_id) {
        return node_id < other.node_id;
    }
    return state_less(state, other.state);
}

bool NodeStateKey::operator==(const NodeStateKey& other) const {
    return node_id == other.node_id && same_state(state, other.state);
}

std::vector<RawPStepPath> enumerate_feasible_raw_psteps(
    const MultiDiGraph& graph,
    int p,
    double time_limit
) {
    require_condition(p >= 1, "p must be at least 1.");
    require_condition(double_greater_or_equal(time_limit, 0.0), "Time limit T must be nonnegative.");

    const std::vector<EdgeRecord>& all_edges = graph.edges();
    const NodeId end_node_id = graph.end_node_id();

    std::vector<RawPStepPath> raw_paths;
    std::vector<int> current_edge_ids;
    std::set<NodeId> visited_physical_nodes;
    int next_path_id = 0;

    const auto dfs = [&](auto&& self,
                         NodeId start_node_id,
                         NodeId current_node_id,
                         State current_state,
                         double current_time) -> void {
        const int current_length = static_cast<int>(current_edge_ids.size());
        if (current_length == 0) {
            return;
        }

        // 길이 조건을 만족하는 중간 경로는 모두 raw p-step 후보로 저장한다.
        if (should_record_raw_path(start_node_id, current_node_id, current_length, p, end_node_id)) {
            raw_paths.push_back(make_raw_path(next_path_id++, current_edge_ids, all_edges));
        }

        if (current_length >= p || current_node_id == end_node_id) {
            return;
        }

        const std::vector<std::size_t>& outgoing_indices =
            graph.outgoing_edge_indices(current_node_id);
        if (outgoing_indices.empty()) {
            return;
        }

        for (std::size_t edge_index : outgoing_indices) {
            const int edge_id = static_cast<int>(edge_index);
            const EdgeRecord& edge = all_edges[edge_index];
            if (!same_state(edge.data.start_state, current_state)) {
                continue;
            }

            const double next_time = current_time + edge.data.time;
            if (!double_less_or_equal(next_time, time_limit)) {
                continue;
            }

            const NodeId next_node_id = edge.v;
            const bool adds_physical_node = graph.is_physical_service_node(next_node_id);
            // 같은 physical service node를 한 raw p-step 안에서 다시 방문하지 않도록 막는다.
            if (adds_physical_node &&
                visited_physical_nodes.find(next_node_id) != visited_physical_nodes.end()) {
                continue;
            }

            current_edge_ids.push_back(edge_id);
            if (adds_physical_node) {
                visited_physical_nodes.insert(next_node_id);
            }

            self(
                self,
                start_node_id,
                next_node_id,
                canonicalize_state(edge.data.end_state),
                next_time
            );

            if (adds_physical_node) {
                visited_physical_nodes.erase(next_node_id);
            }
            current_edge_ids.pop_back();
        }
    };

    // 모든 edge를 시작점으로 삼아 feasible raw p-step을 완전탐색한다.
    // 이 방식은 구현이 단순하고, "start가 depot이면 길이 <= p, 아니면 길이 = p"인
    // 모든 raw p-step을 빠뜨리지 않고 생성하기 쉽다는 장점이 있다.
    // 다만 계산 측면에서는 부분 경로를 여러 번 다시 탐색하는 중복이 생긴다.
    // 예를 들어 raw path가 edge sequence [e3, e7, e9]라면,
    // - e3를 시작 edge로 잡았을 때는 [e3, e7, e9]를 탐색하고
    // - e7를 시작 edge로 잡았을 때는 그 suffix인 [e7, e9]를 다시 탐색한다.
    //
    // 즉 "완전히 같은 raw p-step을 중복 저장"하는 것이 아니라,
    // 서로 다른 시작 edge에서 출발하는 DFS들이 뒤쪽의 공통 suffix를 반복 계산하는 구조다.
    // 특히 (current_node_id, current_state, current_length)가 같아지는 경우에도,
    // visited_physical_nodes의 내용까지 동일하다면 이후 탐색 결과가 사실상 같을 수 있다.
    //
    // 추후 개선 방향:
    // 1. 시작 상태 기준 DP / labeling 방식
    //    edge 단위로 매번 DFS를 새로 시작하는 대신,
    //    동일한 시작 (node, state, length)에서의 확장을 한 번만 계산하는 구조로 바꿀 수 있다.
    // 2. pricing-style 확장 재활용
    //    논문에서의 labeling / shortest path 확장처럼 부분 경로를 label로 관리하면,
    //    동일한 prefix 또는 suffix를 반복 계산하는 비용을 줄일 수 있다.
    for (std::size_t edge_idx = 0; edge_idx < all_edges.size(); ++edge_idx) {
        const EdgeRecord& edge = all_edges[edge_idx];
        if (!double_less_or_equal(edge.data.time, time_limit)) {
            continue;
        }

        current_edge_ids.clear();
        visited_physical_nodes.clear();
        current_edge_ids.push_back(static_cast<int>(edge_idx));

        if (graph.is_physical_service_node(edge.u)) {
            visited_physical_nodes.insert(edge.u);
        }
        if (graph.is_physical_service_node(edge.v)) {
            const auto inserted = visited_physical_nodes.insert(edge.v);
            if (!inserted.second) {
                current_edge_ids.pop_back();
                continue;
            }
        }

        dfs(
            dfs,
            edge.u,
            edge.v,
            canonicalize_state(edge.data.end_state),
            edge.data.time
        );
    }

    return raw_paths;
}

std::vector<CompactPStep> build_compact_psteps(
    const std::vector<RawPStepPath>& raw_paths,
    double time_limit,
    NodeId end_node_id
) {
    std::vector<CompactPStep> compact_psteps;
    int next_pstep_id = 0;

    for (const RawPStepPath& raw_path : raw_paths) {
        // 각 raw path마다 허용되는 극단 tau 값만 사용해 compact p-step을 만든다.
        const std::vector<double> tau_values = extreme_tau_values(raw_path, time_limit, end_node_id);
        for (double tau_value : tau_values) {
            CompactPStep pstep;
            pstep.id = next_pstep_id++;
            pstep.raw_path_id = raw_path.id;
            pstep.edge_ids = raw_path.edge_ids;
            pstep.q = raw_path.q;
            pstep.start_node_id = raw_path.start_node_id;
            pstep.last_node_id = raw_path.last_node_id;
            pstep.start_state = canonicalize_state(raw_path.start_state);
            pstep.last_state = canonicalize_state(raw_path.last_state);
            pstep.total_time = raw_path.total_time;
            pstep.total_cost = raw_path.total_cost;
            pstep.tau = tau_value;
            pstep.node_sequence = raw_path.node_sequence;
            pstep.state_sequence = raw_path.state_sequence;
            compact_psteps.push_back(std::move(pstep));
        }
    }

    return compact_psteps;
}

CompactPStepCoefficients build_compact_pstep_coefficients(
    const MultiDiGraph& graph,
    const std::vector<CompactPStep>& compact_psteps
) {
    const NodeId end_node_id = graph.end_node_id();
    const std::size_t edge_count = graph.number_of_edges();

    CompactPStepCoefficients coefficients;
    coefficients.edge_rows.assign(edge_count, {});
    for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
        coefficients.physical_nodes.push_back(node_id);
    }

    // 각 physical node에서 실제로 등장하는 상태 집합 sigma를 먼저 수집한다.
    std::map<NodeId, std::set<State, decltype(&state_less)>> sigma_sets;
    for (NodeId node_id : coefficients.physical_nodes) {
        sigma_sets.emplace(node_id, std::set<State, decltype(&state_less)>(state_less));
    }

    coefficients.visit_coefficients_by_pstep.reserve(compact_psteps.size());
    coefficients.state_coefficients_by_pstep.reserve(compact_psteps.size());
    coefficients.time_coefficients_by_pstep.reserve(compact_psteps.size());
    coefficients.edge_incidence_by_pstep.reserve(compact_psteps.size());

    for (const CompactPStep& pstep : compact_psteps) {
        if (graph.is_physical_service_node(pstep.start_node_id)) {
            sigma_sets[pstep.start_node_id].insert(canonicalize_state(pstep.start_state));
        }
        if (graph.is_physical_service_node(pstep.last_node_id)) {
            sigma_sets[pstep.last_node_id].insert(canonicalize_state(pstep.last_state));
        }
    }

    for (const auto& entry : sigma_sets) {
        const NodeId node_id = entry.first;
        const std::set<State, decltype(&state_less)>& state_set = entry.second;
        coefficients.sigma_by_node[node_id] = std::vector<State>(state_set.begin(), state_set.end());
    }

    for (const CompactPStep& pstep : compact_psteps) {
        std::map<NodeId, int> visit_coefficients = calculated_visit_coefficients(pstep, graph);
        std::map<NodeStateKey, int> state_coefficients =
            calculated_state_coefficients(pstep, graph);
        std::map<NodeStateKey, double> time_coefficients =
            calculated_time_coefficients(pstep, graph);

        // p-step별 계수와 행 기준 희소 표현을 동시에 구성한다.
        coefficients.visit_coefficients_by_pstep.push_back(visit_coefficients);
        coefficients.state_coefficients_by_pstep.push_back(state_coefficients);
        coefficients.time_coefficients_by_pstep.push_back(time_coefficients);
        coefficients.edge_incidence_by_pstep.push_back(pstep.edge_ids);
        
        //pstep.id는 현재 iteration_num과 동일해서 위에서는 그냥 push_back 하는 것
        for (const auto& entry : visit_coefficients) {
            coefficients.visit_rows[entry.first].push_back({pstep.id, entry.second});
        }
        for (const auto& entry : state_coefficients) {
            coefficients.state_rows[entry.first].push_back({pstep.id, entry.second});
        }
        for (const auto& entry : time_coefficients) {
            coefficients.time_rows[entry.first].push_back({pstep.id, entry.second});
        }
        for (int edge_id : pstep.edge_ids) {
            coefficients.edge_rows[static_cast<std::size_t>(edge_id)].push_back(pstep.id);
        }
    }

    return coefficients;
}

CompactPStepArtifacts build_compact_pstep_artifacts(
    const MultiDiGraph& graph,
    const CompactPStepOptions& options,
    std::ostream* log_stream
) {
    const double time_limit =
        options.time_limit >= 0.0 ? options.time_limit : std::numeric_limits<double>::quiet_NaN();
    require_condition(options.p >= 1, "p must be at least 1.");
    require_condition(
        options.time_limit >= 0.0,
        "Compact p-step artifact builder requires an explicit nonnegative time limit."
    );

    CompactPStepArtifacts artifacts;
    artifacts.p = options.p;
    artifacts.time_limit = time_limit;
    // raw p-step 생성 -> compact p-step 생성 -> 계수 구성 순서로 결과를 만든다.
    artifacts.raw_paths = enumerate_feasible_raw_psteps(graph, options.p, time_limit);

    const NodeId end_node_id = graph.end_node_id();
    artifacts.compact_psteps = build_compact_psteps(artifacts.raw_paths, time_limit, end_node_id);
    artifacts.coefficients =
        build_compact_pstep_coefficients(graph, artifacts.compact_psteps);

    if (options.validate) {
        // 필요 시 생성 결과의 내부 일관성을 즉시 검증한다.
        validate_raw_psteps(graph, artifacts.raw_paths, options.p, time_limit);
        validate_compact_psteps(artifacts.compact_psteps, time_limit, end_node_id);
        validate_compact_pstep_coefficients(
            graph,
            artifacts.compact_psteps,
            artifacts.coefficients
        );
    }

    return artifacts;
}

CompactMasterProblem build_compact_master_problem(
    const MultiDiGraph& graph,
    const CompactPStepArtifacts& artifacts,
    const std::string& log_path
) {
    CompactMasterProblem problem;
    problem.env = std::make_unique<GRBEnv>(true);
    if (!log_path.empty()) {
        problem.env->set(GRB_StringParam_LogFile, log_path);
    }
    problem.env->start();
    problem.model = std::make_unique<GRBModel>(*problem.env);
    problem.model->set(GRB_StringAttr_ModelName, "compact_p_step_master");

    const std::size_t compact_pstep_count = artifacts.compact_psteps.size();
    const std::size_t edge_count = graph.number_of_edges();

    problem.x_vars.reserve(compact_pstep_count);
    // 각 compact p-step에 대해 연속 변수 x_r를 생성한다.
    for (const CompactPStep& pstep : artifacts.compact_psteps) {
        problem.x_vars.push_back(
            problem.model->addVar(
                0.0,
                GRB_INFINITY,
                pstep.total_cost,
                GRB_CONTINUOUS,
                "x_" + std::to_string(pstep.id)
            )
        );
    }

    problem.theta_vars.reserve(edge_count);
    // 각 edge 사용 여부를 나타내는 이진 변수 theta_e를 생성한다.
    for (std::size_t edge_id = 0; edge_id < edge_count; ++edge_id) {
        problem.theta_vars.push_back(
            problem.model->addVar(
                0.0,
                1.0,
                0.0,
                GRB_BINARY,
                "theta_" + std::to_string(edge_id)
            )
        );
    }

    // 각 physical node는 정확히 한 번 들어오고 한 번 나가도록 방문 제약을 둔다.
    for (NodeId node_id : artifacts.coefficients.physical_nodes) {
        GRBLinExpr visit_expr = 0.0;
        const auto found = artifacts.coefficients.visit_rows.find(node_id);
        if (found != artifacts.coefficients.visit_rows.end()) {
            for (const auto& entry : found->second) {
                visit_expr += static_cast<double>(entry.second) * problem.x_vars[entry.first];
            }
        }
        problem.model->addConstr(
            visit_expr == 2.0,
            "visit_" + std::to_string(node_id)
        );
    }

    // 상태 보존 및 시간 관련 제약을 각 (node, sigma)에 대해 추가한다.
    for (const auto& sigma_entry : artifacts.coefficients.sigma_by_node) {
        const NodeId node_id = sigma_entry.first;
        const std::vector<State>& sigma_values = sigma_entry.second;

        for (std::size_t sigma_idx = 0; sigma_idx < sigma_values.size(); ++sigma_idx) {
            const NodeStateKey key{node_id, sigma_values[sigma_idx]};

            GRBLinExpr state_expr = 0.0;
            const auto state_found = artifacts.coefficients.state_rows.find(key);
            if (state_found != artifacts.coefficients.state_rows.end()) {
                for (const auto& entry : state_found->second) {
                    state_expr += static_cast<double>(entry.second) * problem.x_vars[entry.first];
                }
            }
            problem.model->addConstr(
                state_expr == 0.0,
                "state_" + std::to_string(node_id) + "_" + std::to_string(sigma_idx)
            );

            GRBLinExpr time_expr = 0.0;
            const auto time_found = artifacts.coefficients.time_rows.find(key);
            if (time_found != artifacts.coefficients.time_rows.end()) {
                for (const auto& entry : time_found->second) {
                    time_expr += entry.second * problem.x_vars[entry.first];
                }
            }
            problem.model->addConstr(
                time_expr >= 0.0,
                "time_" + std::to_string(node_id) + "_" + std::to_string(sigma_idx)
            );
        }
    }

    // edge와 p-step 선택 사이의 연결 제약을 둔다.
    for (std::size_t edge_id = 0; edge_id < edge_count; ++edge_id) {
        GRBLinExpr edge_expr = 0.0;
        for (int pstep_id : artifacts.coefficients.edge_rows[edge_id]) {
            edge_expr += problem.x_vars[pstep_id];
        }
        problem.model->addConstr(
            edge_expr == problem.theta_vars[edge_id],
            "edge_" + std::to_string(edge_id)
        );
    }

    problem.model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    problem.model->update();

    return problem;
}

void validate_raw_psteps(
    const MultiDiGraph& graph,
    const std::vector<RawPStepPath>& raw_paths,
    int p,
    double time_limit
) {
    const std::vector<EdgeRecord>& all_edges = graph.edges();
    const NodeId end_node_id = graph.end_node_id();

    for (const RawPStepPath& path : raw_paths) {
        // 길이, 시퀀스 크기, 시간 제한 등 기본 구조를 점검한다.
        require_condition(path.q == static_cast<int>(path.edge_ids.size()), "Invalid raw q.");
        require_condition(path.q >= 1, "Raw path must have positive length.");
        require_condition(
            static_cast<int>(path.node_sequence.size()) == path.q + 1,
            "Raw path node sequence length mismatch."
        );
        require_condition(
            static_cast<int>(path.state_sequence.size()) == path.q + 1,
            "Raw path state sequence length mismatch."
        );
        require_condition(
            is_valid_raw_length(path.start_node_id, path.last_node_id, path.q, p, end_node_id),
            "Raw path length condition violated."
        );
        require_condition(
            double_less_or_equal(path.total_time, time_limit),
            "Raw path exceeds time limit T."
        );
        require_condition(path.node_sequence.front() == path.start_node_id, "Raw start node mismatch.");
        require_condition(path.node_sequence.back() == path.last_node_id, "Raw end node mismatch.");
        require_condition(same_state(path.state_sequence.front(), path.start_state), "Raw start state mismatch.");
        require_condition(same_state(path.state_sequence.back(), path.last_state), "Raw end state mismatch.");

        double recomputed_time = 0.0;
        double recomputed_cost = 0.0;
        std::set<NodeId> visited_physical_nodes;

        for (int edge_pos = 0; edge_pos < path.q; ++edge_pos) {
            const int edge_id = path.edge_ids[static_cast<std::size_t>(edge_pos)];
            require_condition(edge_id >= 0, "Raw edge id must be nonnegative.");
            require_condition(
                static_cast<std::size_t>(edge_id) < all_edges.size(),
                "Raw edge id is out of range."
            );

            const EdgeRecord& edge = all_edges[static_cast<std::size_t>(edge_id)];
            recomputed_time += edge.data.time;
            recomputed_cost += edge.data.cost;

            require_condition(
                path.node_sequence[static_cast<std::size_t>(edge_pos)] == edge.u,
                "Raw node sequence is inconsistent with edge source."
            );
            require_condition(
                path.node_sequence[static_cast<std::size_t>(edge_pos + 1)] == edge.v,
                "Raw node sequence is inconsistent with edge target."
            );
            require_condition(
                same_state(path.state_sequence[static_cast<std::size_t>(edge_pos)], edge.data.start_state),
                "Raw state sequence is inconsistent with edge start state."
            );
            require_condition(
                same_state(path.state_sequence[static_cast<std::size_t>(edge_pos + 1)], edge.data.end_state),
                "Raw state sequence is inconsistent with edge end state."
            );

            if (edge_pos + 1 < path.q) {
                // 인접 edge 사이의 노드/상태 연속성을 확인한다.
                const EdgeRecord& next_edge =
                    all_edges[static_cast<std::size_t>(path.edge_ids[static_cast<std::size_t>(edge_pos + 1)])];
                require_condition(edge.v == next_edge.u, "Raw path violates node continuity.");
                require_condition(
                    same_state(edge.data.end_state, next_edge.data.start_state),
                    "Raw path violates state continuity."
                );
            }
        }

        require_condition(double_equal(recomputed_time, path.total_time), "Raw time mismatch.");
        require_condition(double_equal(recomputed_cost, path.total_cost), "Raw cost mismatch.");

        // 하나의 raw p-step 안에서 physical service node 재방문이 없는지 확인한다.
        for (std::size_t pos = 0; pos < path.node_sequence.size(); ++pos) {
            const NodeId node_id = path.node_sequence[pos];
            if (node_id == 0) {
                require_condition(pos == 0U, "Node 0 can appear only at the first position.");
            }
            if (node_id == end_node_id) {
                require_condition(
                    pos + 1U == path.node_sequence.size(),
                    "End node d can appear only at the last position."
                );
            }

            if (graph.is_physical_service_node(node_id)) {
                const auto inserted = visited_physical_nodes.insert(node_id);
                require_condition(
                    inserted.second,
                    "A physical pickup/delivery node appears more than once in a raw p-step."
                );
            }
        }
    }
}

void validate_compact_psteps(
    const std::vector<CompactPStep>& compact_psteps,
    double time_limit,
    NodeId end_node_id
) {
    for (const CompactPStep& pstep : compact_psteps) {
        // compact p-step이 대응 raw path 구조를 그대로 유지하는지 검사한다.
        require_condition(pstep.q >= 1, "Compact p-step must contain at least one edge.");
        require_condition(
            static_cast<int>(pstep.edge_ids.size()) == pstep.q,
            "Compact p-step edge count mismatch."
        );
        require_condition(
            static_cast<int>(pstep.node_sequence.size()) == pstep.q + 1,
            "Compact p-step node sequence length mismatch."
        );
        require_condition(
            static_cast<int>(pstep.state_sequence.size()) == pstep.q + 1,
            "Compact p-step state sequence length mismatch."
        );

        const RawPStepPath raw_path{
            pstep.raw_path_id,
            pstep.edge_ids,
            pstep.q,
            pstep.start_node_id,
            pstep.last_node_id,
            pstep.start_state,
            pstep.last_state,
            pstep.total_time,
            pstep.total_cost,
            pstep.node_sequence,
            pstep.state_sequence,
        };
        const std::vector<double> tau_values = extreme_tau_values(raw_path, time_limit, end_node_id);
        const bool has_valid_tau = std::any_of(
            tau_values.begin(),
            tau_values.end(),
            [&pstep](double tau_value) {
                return double_equal(tau_value, pstep.tau);
            }
        );

        // tau는 해당 raw path에서 허용되는 극단값 중 하나여야 한다.
        require_condition(has_valid_tau, "Compact p-step uses a non-extreme tau value.");
        require_condition(double_greater_or_equal(pstep.tau, 0.0), "tau must be nonnegative.");
        require_condition(
            double_greater_or_equal(pstep.tau + pstep.total_time, 0.0) &&
                double_less_or_equal(pstep.tau + pstep.total_time, time_limit),
            "Compact p-step violates 0 <= tau + t(P) <= T."
        );
    }
}

//같은 방식으로 또 구하는거라 크게 의미 없긴 함
void validate_compact_pstep_coefficients(
    const MultiDiGraph& graph,
    const std::vector<CompactPStep>& compact_psteps,
    const CompactPStepCoefficients& coefficients
) {
    require_condition(
        compact_psteps.size() == coefficients.visit_coefficients_by_pstep.size(),
        "Visit coefficient count does not match the number of compact p-steps."
    );
    require_condition(
        compact_psteps.size() == coefficients.state_coefficients_by_pstep.size(),
        "State coefficient count does not match the number of compact p-steps."
    );
    require_condition(
        compact_psteps.size() == coefficients.time_coefficients_by_pstep.size(),
        "Time coefficient count does not match the number of compact p-steps."
    );

    for (const CompactPStep& pstep : compact_psteps) {
        // 수식에서 정의된 기대 계수와 실제 구성 결과를 비교한다.
        const std::map<NodeId, int> calculated_visit =
            calculated_visit_coefficients(pstep, graph);
        const std::map<NodeStateKey, int> calculated_state =
            calculated_state_coefficients(pstep, graph);
        const std::map<NodeStateKey, double> calculated_time =
            calculated_time_coefficients(pstep, graph);

        const std::map<NodeId, int>& actual_visit =
            coefficients.visit_coefficients_by_pstep[static_cast<std::size_t>(pstep.id)];
        const std::map<NodeStateKey, int>& actual_state =
            coefficients.state_coefficients_by_pstep[static_cast<std::size_t>(pstep.id)];
        const std::map<NodeStateKey, double>& actual_time =
            coefficients.time_coefficients_by_pstep[static_cast<std::size_t>(pstep.id)];

        require_condition(actual_visit == calculated_visit, "a_{i,r} coefficients are inconsistent.");
        require_condition(actual_state == calculated_state, "s_{i,sigma,r} coefficients are inconsistent.");
        require_condition(
            actual_time.size() == calculated_time.size(),
            "q_{i,sigma,r} coefficient count is inconsistent."
        );

        for (const auto& entry : calculated_time) {
            const auto found = actual_time.find(entry.first);
            require_condition(found != actual_time.end(), "Missing q_{i,sigma,r} coefficient.");
            require_condition(
                double_equal(found->second, entry.second),
                "q_{i,sigma,r} coefficient value is inconsistent."
            );
        }

        int positive_count = 0;
        int negative_count = 0;
        for (const auto& entry : actual_state) {
            if (entry.second == 1) {
                ++positive_count;
            } else if (entry.second == -1) {
                ++negative_count;
            } else {
                require_condition(false, "s_{i,sigma,r} must be +1 or -1.");
            }
        }
        require_condition(positive_count <= 1, "s_{i,sigma,r} has more than one +1 entry.");
        require_condition(negative_count <= 1, "s_{i,sigma,r} has more than one -1 entry.");
    }
}

void dump_compact_psteps(
    std::ostream& out,
    const std::vector<CompactPStep>& compact_psteps,
    std::size_t limit
) {
    // 너무 많은 출력이 생기지 않도록 앞부분만 잘라서 보여준다.
    const std::size_t dump_count = std::min(limit, compact_psteps.size());
    out << "[p-step] Dumping " << dump_count << " compact p-steps\n";

    for (std::size_t idx = 0; idx < dump_count; ++idx) {
        const CompactPStep& pstep = compact_psteps[idx];
        out << "  id=" << pstep.id
            << " raw_id=" << pstep.raw_path_id
            << " edges=" << format_edge_sequence(pstep.edge_ids)
            << " nodes=" << format_node_sequence(pstep.node_sequence)
            << " states=" << format_state_sequence(pstep.state_sequence)
            << " tau=" << format_double(pstep.tau)
            << " time=" << format_double(pstep.total_time)
            << " cost=" << format_double(pstep.total_cost)
            << '\n';
    }
}

void solve_compact_master_problem(CompactMasterProblem& problem) {
    require_condition(problem.model != nullptr, "Model must be created before solving.");
    problem.model->optimize();
}

RecoveredSolution recover_incumbent_solution(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const CompactPStepArtifacts& artifacts,
    const CompactMasterProblem& problem
) {
    RecoveredSolution solution;
    solution.runtime_seconds = problem.model->get(GRB_DoubleAttr_Runtime);

    // incumbent가 없으면 복원할 route도 없으므로 빈 해를 그대로 반환한다.
    const int solution_count = problem.model->get(GRB_IntAttr_SolCount);
    if (solution_count <= 0) {
        return solution;
    }

    solution.has_incumbent = true;
    solution.objective_value = problem.model->get(GRB_DoubleAttr_ObjVal);
    solution.active_edge_ids = collect_active_edge_ids(problem);
    const NodeId end_node_id = graph.end_node_id();
    const std::vector<EdgeRecord>& all_edges = graph.edges();
    std::map<NodeId, std::vector<int>> active_outgoing_edges;
    for (int edge_id : solution.active_edge_ids) {
        active_outgoing_edges[all_edges[static_cast<std::size_t>(edge_id)].u].push_back(edge_id);
    }
    for (auto& entry : active_outgoing_edges) {
        std::sort(entry.second.begin(), entry.second.end());
    }

    // active theta support 위에서 아직 어떤 edge가 route에 배정되지 않았는지 추적한다.
    std::vector<bool> is_edge_consumed(all_edges.size(), true);
    for (int edge_id : solution.active_edge_ids) {
        is_edge_consumed[static_cast<std::size_t>(edge_id)] = false;
    }

    // 활성 theta edge들을 실제 feasible route들로 분해해 본다.
    // 이 단계가 성공해야 binary theta 해가 실제 route 집합과 모순 없이 대응된다고 볼 수 있다.
    std::vector<std::vector<int>> route_edge_sets;
    const bool decomposed = decompose_selected_edges(
        data,
        graph,
        all_edges,
        active_outgoing_edges,
        is_edge_consumed,
        route_edge_sets,
        end_node_id
    );
    if (!decomposed) {
        // 분해 실패 시 활성 theta와 양의 x 값을 함께 남겨 디버깅 가능하게 한다.
        std::ostringstream debug;
        debug << "The selected theta-edge set could not be decomposed into feasible routes.";
        debug << " Active theta: ";
        for (std::size_t idx = 0; idx < solution.active_edge_ids.size(); ++idx) {
            if (idx > 0U) {
                debug << " || ";
            }
            const int edge_id = solution.active_edge_ids[idx];
            debug << format_debug_edge(edge_id, all_edges[static_cast<std::size_t>(edge_id)]);
        }
        debug << " Positive x: ";
        bool has_positive_x = false;
        for (std::size_t idx = 0; idx < artifacts.compact_psteps.size(); ++idx) {
            const double value = problem.x_vars[idx].get(GRB_DoubleAttr_X);
            if (value <= 1e-6) {
                continue;
            }
            if (has_positive_x) {
                debug << " || ";
            }
            has_positive_x = true;
            debug << format_debug_pstep(artifacts.compact_psteps[idx], value);
        }
        throw std::runtime_error(debug.str());
    }

    // edge 단위 route를 사람이 읽을 수 있는 route 출력 형식으로 변환한다.
    int next_route_id = 0;
    std::set<int> recovered_edge_ids;
    std::set<NodeId> recovered_service_nodes;
    for (const std::vector<int>& route_edge_ids : route_edge_sets) {
        for (int edge_id : route_edge_ids) {
            recovered_edge_ids.insert(edge_id);

            const NodeId arrival_node_id = all_edges[static_cast<std::size_t>(edge_id)].v;
            if (graph.is_physical_service_node(arrival_node_id)) {
                const auto inserted = recovered_service_nodes.insert(arrival_node_id);
                require_condition(
                    inserted.second,
                    "Recovered solution visits the same service node more than once."
                );
            }
        }

        RecoveredRouteSolution route = build_route_solution(next_route_id, route_edge_ids, data, graph);
        if (!route.actions.empty()) {
            route.route_id = next_route_id++;
            solution.routes.push_back(std::move(route));
        }
    }

    // 복원된 route들이 활성 theta support와 정확히 일치하는지 마지막으로 확인한다.
    // 일부 edge가 누락되거나 중복 사용되면 복원 실패로 간주한다.
    for (int edge_id : solution.active_edge_ids) {
        require_condition(
            is_edge_consumed[static_cast<std::size_t>(edge_id)],
            "The selected theta-edge set could not be fully decomposed into routes."
        );
    }

    require_condition(
        recovered_edge_ids.size() == solution.active_edge_ids.size(),
        "Recovered route edges do not match the active theta support size."
    );
    for (int edge_id : solution.active_edge_ids) {
        require_condition(
            recovered_edge_ids.find(edge_id) != recovered_edge_ids.end(),
            "Recovered route edges do not match the active theta support."
        );
    }
    for (NodeId node_id = 1; node_id < end_node_id; ++node_id) {
        require_condition(
            recovered_service_nodes.find(node_id) != recovered_service_nodes.end(),
            "Recovered solution does not cover every service node."
        );
    }

    return solution;
}

void write_recovered_solution(std::ostream& out, const RecoveredSolution& solution) {
    out << "Solution:\n";

    if (!solution.has_incumbent) {
        out << "  No incumbent solution available\n";
        out << "Solution done \n";
        return;
    }

    // route별 비용은 이미 복원되어 있으므로 전체 비용은 route cost 합으로 출력한다.
    double total_route_cost = 0.0;
    for (const RecoveredRouteSolution& route : solution.routes) {
        total_route_cost += route.total_cost;
    }

    out << "  Number of routes " << solution.routes.size() << '\n';
    out << "  Total cost " << format_solution_value(total_route_cost) << '\n';

    // 각 route를 Path, Loc, Patt, Load 형식으로 정리해 출력한다.
    for (const RecoveredRouteSolution& route : solution.routes) {
        out << "Route " << route.route_id << ":\n";
        out << "  Cost " << format_solution_value(route.total_cost) << '\n';
        out << "  Time " << format_solution_value(route.total_time) << '\n';

        out << "  Path :";
        for (const RouteSolutionAction& action : route.actions) {
            out << ' ' << action.request_id;
        }
        out << " \n";

        out << "  Loc  :";
        for (const RouteSolutionAction& action : route.actions) {
            out << ' ' << action.location;
        }
        out << " \n";

        out << "  Patt :";
        for (const RouteSolutionAction& action : route.actions) {
            out << ' ' << action.pattern;
        }
        out << " \n";

        out << "  Load :";
        for (const RouteSolutionAction& action : route.actions) {
            out << ' ' << action.load;
        }
        out << " \n";
    }

    out << "Solution done \n";
}

}  // namespace spdp
