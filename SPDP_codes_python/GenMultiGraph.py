from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple
import argparse

import networkx as nx

try:
    from ReadData import SPDPData, read_spdp_data
except ImportError:  # pragma: no cover
    from .ReadData import SPDPData, read_spdp_data


NodeId = int
StateToken = Tuple[str, int, int]  # ('N',-1,-1), ('E', h, -1), ('F', h, l)
State = Tuple[StateToken, StateToken]  # capacity = 2 (multiset encoded as sorted tuple)

N_TOKEN: StateToken = ("N", -1, -1)


@dataclass(frozen=True)
class NodeSpec:
    node_id: NodeId
    kind: str  # start, end, pickup, delivery
    location: int
    request_idx: Optional[int] = None
    container_type: Optional[int] = None
    landfill_location: Optional[int] = None


def _make_e(container_type: int) -> StateToken:
    return ("E", container_type, -1)


def _make_f(container_type: int, landfill_location: int) -> StateToken:
    return ("F", container_type, landfill_location)


def _canonical_state(tokens: Iterable[StateToken]) -> State: ##state 순서를 일관되게 하여 같은 state지만 순서만 다른 state가 나오지 않게 함
    ordered = tuple(sorted(tokens))
    if len(ordered) != 2:
        raise ValueError(f"State must have exactly 2 slots, got {len(ordered)}")
    return ordered  # type: ignore[return-value]


def _state_to_str(state: State) -> str:
    return f"{state[0]}|{state[1]}"


def _build_node_specs(data: SPDPData) -> Dict[NodeId, NodeSpec]:
    n_requests = len(data.requests)
    end_depot_id = 2 * n_requests + 1
    virtual_location = int(data.locations)
    nodes: Dict[NodeId, NodeSpec] = {
        0: NodeSpec(node_id=0, kind="start", location=virtual_location),
        end_depot_id: NodeSpec(node_id=end_depot_id, kind="end", location=virtual_location),
    }
    for idx, req in enumerate(data.requests):
        pickup_id = idx + 1
        delivery_id = n_requests + idx + 1
        nodes[pickup_id] = NodeSpec(
            node_id=pickup_id,
            kind="pickup",
            location=req.from_id,
            request_idx=idx,
            container_type=req.container_type,
            landfill_location=req.to_id,
        )
        nodes[delivery_id] = NodeSpec(
            node_id=delivery_id,
            kind="delivery",
            location=req.from_id,
            request_idx=idx,
            container_type=req.container_type,
            landfill_location=req.to_id,
        )
    return nodes


def _generate_state_space(data: SPDPData) -> List[State]:
    container_types = sorted({req.container_type for req in data.requests})
    full_items = sorted({(req.container_type, req.to_id) for req in data.requests})

    atom_space: List[StateToken] = [N_TOKEN]
    atom_space.extend(_make_e(h) for h in container_types)
    atom_space.extend(_make_f(h, l) for h, l in full_items)

    states: List[State] = []
    for pair in combinations_with_replacement(atom_space, 2):
        states.append(_canonical_state(pair))
    return states


def _service_time(node: NodeSpec, data: SPDPData) -> float:
    if node.kind == "pickup":
        return float(data.time_pickup)
    if node.kind == "delivery":
        return float(data.time_delivery)
    return 0.0


def _apply_service(node: NodeSpec, state: State) -> Optional[State]:
    tokens = [state[0], state[1]]
    empty_state = _canonical_state((N_TOKEN, N_TOKEN))

    if node.kind == "start":
        return _canonical_state(tokens)
    if node.kind == "end":
        return _canonical_state(tokens) if _canonical_state(tokens) == empty_state else None

    if node.kind == "pickup":
        for i, token in enumerate(tokens):
            if token[0] == "N":
                tokens[i] = _make_f(node.container_type, node.landfill_location)  # type: ignore[arg-type]
                return _canonical_state(tokens)
        return None

    if node.kind == "delivery":
        target = _make_e(node.container_type)  # type: ignore[arg-type]
        for i, token in enumerate(tokens):
            if token == target:
                tokens[i] = N_TOKEN
                return _canonical_state(tokens)
        return None

    raise ValueError(f"Unknown node kind: {node.kind}")


def _candidate_sequences(state: State) -> List[Tuple[int, ...]]:
    full_locations = [token[2] for token in state if token[0] == "F"]
    if len(full_locations) == 0:
        return [tuple()]
    if len(full_locations) == 1:
        loc = full_locations[0]
        return [tuple(), (loc,)]
    if len(full_locations) == 2:
        l1, l2 = full_locations
        if l1 == l2:
            return [tuple(), (l1,)]
        candidates = [tuple(), (l1,), (l2,), (l1, l2), (l2, l1)]
        unique: List[Tuple[int, ...]] = []
        seen: Set[Tuple[int, ...]] = set()
        for seq in candidates:
            if seq not in seen:
                seen.add(seq)
                unique.append(seq)
        return unique
    raise ValueError("Capacity is 2, so full skip count must be <= 2")


def _apply_emptying(state: State, sequence_pi: Sequence[int]) -> Tuple[State, int]:
    tokens = [state[0], state[1]]
    emptied_count = 0
    for treatment_location in sequence_pi:
        for i, token in enumerate(tokens):
            if token[0] == "F" and token[2] == treatment_location:
                tokens[i] = _make_e(token[1])
                emptied_count += 1
    return _canonical_state(tokens), emptied_count


def _path_metric(metric_matrix, start_loc: int, sequence_pi: Sequence[int], end_loc: int) -> float:
    points = [start_loc, *sequence_pi, end_loc]
    value = 0.0
    for i in range(len(points) - 1):
        value += float(metric_matrix[points[i], points[i + 1]])
    return value


def _dominates(a: Dict, b: Dict) -> bool:
    return (a["time"] <= b["time"] and a["cost"] <= b["cost"]) and (
        a["time"] < b["time"] or a["cost"] < b["cost"]
    )


def _violates_single_request_state_rules(
    state: State,
    singleton_types: Set[int],
    singleton_full_pairs: Set[Tuple[int, int]],
) -> bool:
    # Rule 1: if a container type appears in only one request, states containing two of that type are forbidden.
    type_count: Dict[int, int] = {}
    for token in state:
        if token[0] in {"E", "F"}:
            type_count[token[1]] = type_count.get(token[1], 0) + 1
    for container_type, count in type_count.items():
        if count >= 2 and container_type in singleton_types:
            return True

    # Rule 2: if a (type, location) full request appears only once, duplicated full-full state is forbidden.
    if state[0][0] == "F" and state[1][0] == "F":
        if state[0][1] == state[1][1] and state[0][2] == state[1][2]:
            if (state[0][1], state[0][2]) in singleton_full_pairs:
                return True

    return False


def build_multigraph(
    data: SPDPData,
    prune_dominated_edges: bool = True,
    logger: Callable[[str], None] = print,
) -> nx.MultiDiGraph:
    graph = nx.MultiDiGraph()
    node_specs = _build_node_specs(data)
    virtual_location = int(data.locations)
    all_states = _generate_state_space(data)
    empty_state = _canonical_state((N_TOKEN, N_TOKEN))
    request_type_count: Dict[int, int] = {}
    request_full_pair_count: Dict[Tuple[int, int], int] = {}
    for req in data.requests:
        request_type_count[req.container_type] = request_type_count.get(req.container_type, 0) + 1
        pair_key = (req.container_type, req.to_id)
        request_full_pair_count[pair_key] = request_full_pair_count.get(pair_key, 0) + 1
    singleton_types = {container_type for container_type, count in request_type_count.items() if count == 1}
    singleton_full_pairs = {pair for pair, count in request_full_pair_count.items() if count == 1}

    for node in node_specs.values():
        graph.add_node(
            node.node_id,
            kind=node.kind,
            location=node.location,
            request_idx=node.request_idx,
            container_type=node.container_type,
            landfill_location=node.landfill_location,
        )

    feasible_after_states: Dict[NodeId, Set[State]] = {}
    for node in node_specs.values():
        states = set()
        for state in all_states:
            after = _apply_service(node, state)
            if after is not None:
                states.add(after)
        feasible_after_states[node.node_id] = states

    # Fixed rule: virtual start (0) must be empty state.
    feasible_after_states[0] = {empty_state}

    edge_bucket: Dict[Tuple[NodeId, NodeId, State, State], List[Dict]] = {}
    raw_edge_count = 0

    #u->v edge
    for u in node_specs.values():
        if u.kind == "end":
            continue
        for v in node_specs.values():
            if v.kind == "start":
                continue
            if u.node_id == v.node_id:
                continue

            for sigma_u in feasible_after_states[u.node_id]:
                for sequence_pi in _candidate_sequences(sigma_u):
                    sigma_arrival, emptied_count = _apply_emptying(sigma_u, sequence_pi)
                    sigma_v = _apply_service(v, sigma_arrival)
                    if sigma_v is None:
                        continue
                    if _violates_single_request_state_rules(
                        sigma_u, singleton_types, singleton_full_pairs
                    ) or _violates_single_request_state_rules(
                        sigma_v, singleton_types, singleton_full_pairs
                    ):
                        continue

                    travel_time = _path_metric(data.time, u.location, sequence_pi, v.location)
                    travel_cost = _path_metric(data.distance, u.location, sequence_pi, v.location)
                    if u.node_id == 0:
                        travel_cost += float(data.fixed_vehicle_cost)
                    total_time = (
                        travel_time
                        + emptied_count * float(data.time_empty)
                        + _service_time(v, data)
                    )

                    if total_time > float(data.time_limit):
                        continue
                    
                    #u서비스 이후 -> v서비스 이후를 의미
                    edge_data = {
                        "sequence_pi": tuple(sequence_pi),
                        "time": float(total_time),
                        "cost": float(travel_cost),
                        "start_state": sigma_u,
                        "end_state": sigma_v,
                        "start_state_str": _state_to_str(sigma_u),
                        "end_state_str": _state_to_str(sigma_v),
                    }
                    key = (u.node_id, v.node_id, sigma_u, sigma_v)
                    edge_bucket.setdefault(key, []).append(edge_data)
                    raw_edge_count += 1

    kept_edge_count = 0
    for (u_id, v_id, _, _), candidates in edge_bucket.items():
        if prune_dominated_edges:
            kept: List[Dict] = []
            for i, candidate in enumerate(candidates):
                dominated = False
                for j, other in enumerate(candidates):
                    if i == j:
                        continue
                    if _dominates(other, candidate):
                        dominated = True
                        break
                if not dominated:
                    kept.append(candidate)
        else:
            kept = candidates

        for edge_data in kept:
            graph.add_edge(u_id, v_id, **edge_data)
            kept_edge_count += 1

    logger(f"[GenMultiGraph] Node count: {graph.number_of_nodes()}")
    logger(f"[GenMultiGraph] State count: {len(all_states)}")
    logger(f"[GenMultiGraph] Matrix shapes: time={data.time.shape}, distance={data.distance.shape}")
    logger(f"[GenMultiGraph] Virtual location index: {virtual_location}")
    logger(f"[GenMultiGraph] Raw edge count: {raw_edge_count}")
    if prune_dominated_edges:
        logger(f"[GenMultiGraph] Pareto-pruned edge count: {kept_edge_count}")
    else:
        logger(f"[GenMultiGraph] Pruning disabled edge count: {kept_edge_count}")
    logger(f"[GenMultiGraph] Graph edge count: {graph.number_of_edges()}")

    return graph


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate internalized treatment multi-graph.")
    parser.add_argument(
        "instance",
        nargs="?",
        default="RecDep_day_A1.dat",
        help="Instance filename (e.g., RecDep_day_A1.dat) or full path",
    )
    args = parser.parse_args()

    data = read_spdp_data(args.instance)
    build_multigraph(data)
