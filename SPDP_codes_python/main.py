from __future__ import annotations

import argparse

try:
    from ReadData import read_spdp_data
    from GenMultiGraph import build_multigraph
except ImportError:  # pragma: no cover
    from .ReadData import read_spdp_data
    from .GenMultiGraph import build_multigraph


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Load SPDP instance and build Internalized Treatment Multi-graph."
    )
    parser.add_argument(
        "instance",
        nargs="?",
        default="A0.dat",
        help="Instance filename (e.g., RecDep_day_A1.dat) or full path",
    )
    parser.add_argument(
        "--sample-edges",
        type=int,
        default=5,
        help="Number of sample edges to print.",
    )
    parser.add_argument(
        "--prune-dominated-edges",
        type=int,
        choices=[0, 1],
        default=0,
        help="1: enable Pareto dominated-edge pruning, 0: disable pruning.",
    )
    args = parser.parse_args()

    data = read_spdp_data(args.instance)
    graph = build_multigraph(
        data,
        prune_dominated_edges=bool(args.prune_dominated_edges),
    )

    print(f"[main] Loaded instance: {args.instance}")
    print(f"[main] Requests: {len(data.requests)}")
    print(f"[main] Locations: {data.locations}")
    print(f"[main] Node count: {graph.number_of_nodes()}")
    print(f"[main] Edge count: {graph.number_of_edges()}")

    if args.sample_edges > 0:
        print(f"[main] Sample edges (up to {args.sample_edges}):")
        shown = 0
        for u, v, edge_data in graph.edges(data=True):
            print(
                f"  {u} -> {v} | pi={edge_data['sequence_pi']} "
                f"| time={edge_data['time']:.2f} | cost={edge_data['cost']:.2f}"
            )
            shown += 1
            if shown >= args.sample_edges:
                break

    # Validation print for A0-style inspection: print all edges on selected node pairs.
    target_pairs = [(0, 1),(1, 2), (1, 3), (3, 1), (3, 4)]
    print("[main] Full edge list for selected node pairs:")
    for u, v in target_pairs:
        pair_edges = list(graph.edges(nbunch=[u], keys=True, data=True))
        pair_edges = [e for e in pair_edges if e[1] == v]
        print(f"  Pair ({u}, {v}) -> {len(pair_edges)} edges")
        for _, _, edge_key, edge_data in pair_edges:
            print(
                f"    key={edge_key} | pi={edge_data['sequence_pi']} "
                f"| time={edge_data['time']:.2f} | cost={edge_data['cost']:.2f} "
                f"| start={edge_data['start_state_str']} | end={edge_data['end_state_str']}"
            )


if __name__ == "__main__":
    main()
