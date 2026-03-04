from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Union
import re

import numpy as np


Number = Union[int, float]


@dataclass
class Request:
    # One transportation request r:
    # - from_id: pickup location index alpha(r) (0-based, matrix index와 직접 매칭)
    # - to_id: treatment location index beta(r) (0-based, matrix index와 직접 매칭)
    # - container_type: container/skip type phi(r)
    from_id: int
    to_id: int
    container_type: int


@dataclass
class SPDPData:
    # Parsed SPDP instance data structure used for research code:
    # 1) global parameters
    # 2) requests list
    # 3) distance/time matrices
    #
    # Global parameters
    # - fixed_vehicle_cost: fixed cost per vehicle
    # - time_pickup: pickup service time s^p
    # - time_empty: emptying service time s^e
    # - time_delivery: delivery service time s^d
    # - time_limit: route time limit T per vehicle
    # - locations: number of physical locations (matrix size)
    #
    # Request block
    # - requests: List[Request]
    #
    # Matrix block
    # - distance[i, j]: travel cost c_ij
    # - time[i, j]: travel time t_ij
    # Both are numpy arrays of shape (locations, locations), 0-based indices.
    fixed_vehicle_cost: Number
    time_pickup: Number
    time_empty: Number
    time_delivery: Number
    time_limit: Number
    locations: int
    requests: List[Request]
    distance: np.ndarray
    time: np.ndarray


def _split_tokens(line: str) -> List[str]:
    return re.split(r"\s+", line.strip())


def _parse_number(token: str) -> Number:
    value = float(token)
    return int(value) if value.is_integer() else value


def _resolve_data_path(filename: str) -> Path:
    input_path = Path(filename)
    if input_path.is_file():
        return input_path.resolve()

    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    candidates = [
        project_root / "SPDP_data" / filename,
        Path.cwd() / "SPDP_data" / filename,
        Path.cwd() / filename,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(f"Could not find data file: {filename}")


def _parse_matrix(lines: List[str], start_idx: int, locations: int, section_name: str) -> tuple[np.ndarray, int]:
    rows = []
    idx = start_idx
    while idx < len(lines) and len(rows) < locations:
        line = lines[idx].strip()
        if not line:
            idx += 1
            continue
        tokens = _split_tokens(line)
        try:
            row = [float(token) for token in tokens]
        except ValueError as exc:
            raise ValueError(
                f"Non-numeric token found while reading {section_name} matrix at line {idx + 1}: {line}"
            ) from exc

        if len(row) != locations:
            raise ValueError(
                f"{section_name} matrix row length mismatch at line {idx + 1}: "
                f"expected {locations}, got {len(row)}"
            )
        rows.append(row)
        idx += 1

    if len(rows) != locations:
        raise ValueError(
            f"{section_name} matrix row count mismatch: expected {locations}, got {len(rows)}"
        )

    return np.array(rows, dtype=float), idx


def read_spdp_data(filename: str) -> SPDPData:
    # Return one fully parsed SPDPData object from a .dat filename.
    data_path = _resolve_data_path(filename)
    with data_path.open("r", encoding="utf-8") as f:
        raw_lines = f.readlines()

    params = {
        "FixedVehicleCost": None,
        "TimePickUp": None,
        "TimeEmpty": None,
        "TimeDelivery": None,
        "TimeLimit": None,
        "LOCATIONS": None,
    }

    requests: List[Request] = []
    request_count = None
    distance = None
    time = None

    idx = 0
    section = "header"
    while idx < len(raw_lines):
        line = raw_lines[idx].strip()
        if not line:
            idx += 1
            continue

        tokens = _split_tokens(line)
        key = tokens[0]
        key_lower = key.lower()

        if key_lower == "requests":
            section = "requests"
            if len(tokens) > 1:
                request_count = int(float(tokens[1]))
            idx += 1
            continue

        if key_lower == "distance":
            if params["LOCATIONS"] is None:
                raise ValueError("LOCATIONS must be parsed before Distance matrix.")
            distance, idx = _parse_matrix(raw_lines, idx + 1, int(params["LOCATIONS"]), "Distance")
            section = "distance"
            continue

        if key_lower == "time":
            if params["LOCATIONS"] is None:
                raise ValueError("LOCATIONS must be parsed before Time matrix.")
            time, idx = _parse_matrix(raw_lines, idx + 1, int(params["LOCATIONS"]), "Time")
            section = "time"
            continue

        if section == "header":
            if key in params and len(tokens) > 1:
                params[key] = _parse_number(tokens[1])
        elif section == "requests":
            if key.upper() == "ID":
                idx += 1
                continue
            if len(tokens) < 6:
                raise ValueError(f"Invalid request row format at line {idx + 1}: {line}")

            from_id = int(float(tokens[2]))
            container_type = int(float(tokens[3]))
            to_id = int(float(tokens[5]))
            requests.append(
                Request(from_id=from_id, to_id=to_id, container_type=container_type)
            )

        idx += 1

    missing_params = [k for k, v in params.items() if v is None]
    if missing_params:
        raise ValueError(f"Missing required global parameters: {missing_params}")
    if distance is None:
        raise ValueError("Distance matrix was not parsed.")
    if time is None:
        raise ValueError("Time matrix was not parsed.")
    if request_count is not None and request_count != len(requests):
        raise ValueError(
            f"REQUESTS count mismatch: header={request_count}, parsed={len(requests)}"
        )

    return SPDPData(
        fixed_vehicle_cost=params["FixedVehicleCost"],
        time_pickup=params["TimePickUp"],
        time_empty=params["TimeEmpty"],
        time_delivery=params["TimeDelivery"],
        time_limit=params["TimeLimit"],
        locations=int(params["LOCATIONS"]),
        requests=requests,
        distance=distance,
        time=time,
    )


if __name__ == "__main__":
    sample_file = "RecDep_day_A1.dat"
    parsed = read_spdp_data(sample_file)

    print(f"Loaded file: {sample_file}")
    print("Global Parameters:")
    print(f"  FixedVehicleCost: {parsed.fixed_vehicle_cost}")
    print(f"  TimePickUp: {parsed.time_pickup}")
    print(f"  TimeEmpty: {parsed.time_empty}")
    print(f"  TimeDelivery: {parsed.time_delivery}")
    print(f"  TimeLimit: {parsed.time_limit}")
    print(f"  LOCATIONS: {parsed.locations}")

    if parsed.requests:
        first = parsed.requests[0]
        print("First Request:")
        print(f"  FROM_ID: {first.from_id}")
        print(f"  TO_ID: {first.to_id}")
        print(f"  CONTAINER_TYPE: {first.container_type}")
    else:
        print("First Request: None")

    print(f"Distance matrix shape: {parsed.distance.shape}")
    print(f"Time matrix shape: {parsed.time.shape}")
