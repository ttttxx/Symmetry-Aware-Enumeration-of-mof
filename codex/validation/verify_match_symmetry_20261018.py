"""Independent finite-group checks for the MATCH appendix (no coloring DFS).

Run with .venv/bin/python codex/validation/verify_match_symmetry_20261018.py.
Uses explicit finite subgroup closure and integer fixed-set subtraction, not
production subgroup generation, canonical pruning, or quota enumeration.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from main import load_graph, enrich_graph_data_with_topology
from symmetry import SymmetryEngine


def certify():
    graph = enrich_graph_data_with_topology(load_graph(str(ROOT / 'data/vba/vba_graph.json')))
    perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
    size = len(perms)
    n = graph['n_nodes']
    assert size == 16 and n == 32
    lookup = {p: i for i, p in enumerate(perms)}
    identity = lookup[tuple(range(n))]
    mult = [[lookup[tuple(p[q[v]] for v in range(n))] for q in perms] for p in perms]
    degrees = [node['degree'] for node in graph['nodes']]
    assert set(degrees) == {4}
    assert all(degrees[v] == degrees[p[v]] for p in perms for v in range(n))

    def closure(generators):
        current = set(generators) | {identity}
        while True:
            expanded = current | {mult[a][b] for a in current for b in current}
            if expanded == current:
                return frozenset(current)
            current = expanded

    # Close each subgroup with every possible group element until stable.
    subgroups = [frozenset({identity})]
    seen = set(subgroups)
    for subgroup in subgroups:
        for element in range(size):
            joined = closure(set(subgroup) | {element})
            if joined not in seen:
                seen.add(joined)
                subgroups.append(joined)
    subgroups.sort(key=lambda subgroup: (-len(subgroup), tuple(sorted(subgroup))))

    def partition(subgroup):
        remaining = set(range(n))
        blocks = []
        while remaining:
            v = min(remaining)
            block = frozenset(perms[g][v] for g in subgroup)
            remaining -= block
            blocks.append(tuple(sorted(block)))
        return tuple(blocks)

    partitions = {h: partition(h) for h in subgroups}
    expected_tasks = {1: 1, 2: 15, 4: 35, 8: 15, 16: 1}
    assert dict(Counter(map(len, subgroups))) == expected_tasks
    for order in expected_tasks:
        assert len({partitions[h] for h in subgroups if len(h) == order}) == expected_tasks[order]
    # These benchmarks have two regular vertex orbits of a C2^4 action.
    assert all(mult[g][g] == identity for g in range(size))
    assert all(mult[g][h] == mult[h][g] for g in range(size) for h in range(size))
    assert all(len(block) == len(h) for h in subgroups for block in partitions[h])

    certificates = {}
    for q, pool in [(2, 'C4x2_Cu_B'), (3, 'C4x3_Cu_B_Cd')]:
        # Fix(H) is partitioned by its exact stabilizer supergroups.
        exact = {}
        for h in subgroups:
            exact[h] = q ** len(partitions[h]) - sum(exact[k] for k in exact if h < k)
            assert exact[h] >= 0
        assert sum(exact.values()) == q ** n
        layers = Counter()
        for h, count in exact.items():
            layers[len(h)] += count
        counts = {}
        for order, labeled_count in layers.items():
            assert labeled_count * order % size == 0
            counts[order] = labeled_count * order // size
        # Independent Burnside count from individual element cycles.
        burnside = sum(q ** len(partition(closure({g}))) for g in range(size))
        assert burnside % size == 0 and sum(counts.values()) == burnside // size
        assert sum(Fraction(count, order) for order, count in counts.items()) == Fraction(q ** n, size)
        record = json.loads((ROOT / f'results/vba_rho_stratification_20261010/{pool}_rho_quota_summary.json').read_text())
        assert record['output_cap'] >= 5 * record['per_level_quota']
        for layer in record['layer_report']:
            order = layer['stabilizer_order']
            assert layer['task_timeouts'] == 0
            assert layer['fixed_layer_tasks'] == expected_tasks[order]
            assert layer['returned'] == min(record['per_level_quota'], counts[order])
            if layer['status'] == 'layer_exhausted':
                assert layer['returned'] == counts[order]
        certificates[pool] = {
            'exact_orbits_by_stabilizer_order': dict(sorted(counts.items(), reverse=True)),
            'total_orbits': sum(counts.values()),
            'quota_contract_matches_archived_run': True,
        }
    return {'group_order': size, 'vertices': n, 'subgroups_by_order': expected_tasks,
            'checks': certificates}


if __name__ == '__main__':
    result = certify()
    destination = Path(__file__).with_suffix('.json')
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
