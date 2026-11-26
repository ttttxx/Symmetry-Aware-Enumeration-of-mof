"""Inputs and independent counting checks for the 2026-11-16 scale study."""
import collections
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASES = [("ffc",2,2),("ffc",3,3),("ffc",4,4),("ffc",5,5),
         ("ffc",5,8),("ffc",5,12),("ffc",5,15),
         ("fjh",1,2),("fjh",2,1),("fjh",1,3),("fjh",1,4),
         ("lil",2,2),("lil",3,3),("lil",3,4),("lil",4,3)]


def pools_for(c3, c4):
    pools = {}
    for degree, count in [(3,c3),(4,c4)]:
        available = sorted(path.stem for path in (ROOT/"nodes").glob(f"{degree}c_*.cif"))
        if not 1 <= count <= len(available):
            raise ValueError(f"Requested {count} labels for degree {degree}, available {len(available)}")
        pools[degree] = available[:count]
    return pools


def orbit_count(permutations, degrees, pools):
    total = 0
    for permutation in permutations:
        visited = set()
        fixed = 1
        for start in range(len(degrees)):
            if start in visited:
                continue
            current = start
            while current not in visited:
                if degrees[current] != degrees[start]:
                    raise ValueError("Group operation mixes degree-dependent label pools")
                visited.add(current)
                current = permutation[current]
            fixed *= len(pools[degrees[start]])
        total += fixed
    result, remainder = divmod(total, len(permutations))
    assert remainder == 0
    return result


def prepare(topology, c3, c4):
    import numpy as np
    from pymatgen.core import Structure
    import libcasm.xtal as xtal
    import libcasm.configuration as config
    from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
    from symmetry import SymmetryEngine

    pools = pools_for(c3,c4)
    graph = enrich_graph_data_with_topology(load_graph(str(ROOT/"data"/topology/f"{topology}_graph.json")))
    degrees = [node["degree"] for node in graph["nodes"]]
    perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
    inverse_set = {tuple(p.index(i) for i in range(len(p))) for p in perms}
    structure = Structure.from_file(ROOT/"data"/"mof"/f"{topology}.cif")
    assert len(structure) == len(degrees)
    assert all(np.linalg.norm(site.frac_coords-np.array(node["frac"])) < 1e-3
               for site,node in zip(structure,graph["nodes"]))
    labels_by_degree = {degree:[f"L{degree}_{i:03}" for i in range(len(pool))]
                        for degree,pool in pools.items()}
    occ = [labels_by_degree[degree] for degree in degrees]
    occupants = {label:xtal.Occupant(label) for labels in labels_by_degree.values() for label in labels}
    prim = config.Prim(xtal.Prim(xtal.Lattice(structure.lattice.matrix),
                       np.array([site.frac_coords for site in structure]).T, occ, occupants=occupants))
    cell = config.Supercell(prim, transformation_matrix_to_super=np.eye(3,dtype=int))
    proper, full = [], []
    op,end = config.SupercellSymOp.begin(cell),config.SupercellSymOp.end(cell)
    while op != end:
        copied = op.copy()
        full.append(copied)
        if np.linalg.det(op.to_symop().matrix()) > 0:
            proper.append(copied)
        op.next()
    proper_perms = [tuple(op.permute_index(i) for i in range(len(degrees))) for op in proper]
    full_perms = [tuple(op.permute_index(i) for i in range(len(degrees))) for op in full]
    if set(proper_perms) != inverse_set or len(proper_perms) != len(perms):
        raise ValueError(f"{topology}: proper CASM site action differs from ours")
    naive = math.prod(len(pools[degree]) for degree in degrees)
    meta = dict(topology=topology,c3=c3,c4=c4,sites=len(degrees),
                degrees=dict(collections.Counter(degrees)),group_size=len(perms),
                casm_full_size=len(full),naive=naive,
                expected_orbits=orbit_count(perms,degrees,pools),
                expected_full_orbits=orbit_count(full_perms,degrees,pools),
                pools=pools,site_action_matches=True)
    return dict(graph=graph,degrees=degrees,perms=perms,pools=pools,
                requirements=calculate_topology_requirements(graph),
                prim=prim,cell=cell,proper=proper,meta=meta)


def record_digest(records, labels, perms):
    """Normalize CASM max/ours min outside timing, then hash all sorted records."""
    import hashlib
    import numpy as np
    codes = {label:i for i,label in enumerate(labels)}
    n = len(perms[0])
    out = np.empty((len(records),n+8),dtype=np.uint8)
    inverses = [tuple(p.index(i) for i in range(n)) for p in perms]
    for start in range(0,len(records),8192):
        chunk = records[start:start+8192]
        data = np.array([[codes[label] for label in row["coloring"]] for row in chunk],dtype=np.uint8)
        best = data.copy()
        indices = np.arange(len(chunk))
        for inverse in inverses:
            image = data[:,inverse]
            differences = image != best
            first = differences.argmax(axis=1)
            smaller = differences.any(axis=1) & (image[indices,first] < best[indices,first])
            best[smaller] = image[smaller]
        stop = start+len(chunk)
        out[start:stop,:n] = best
        for column,key in [(n,"stabilizer_size"),(n+4,"orbit_size")]:
            values = np.array([row[key] for row in chunk],dtype=">u4")
            out[start:stop,column:column+4] = values.view(np.uint8).reshape(-1,4)
    view = out.view(np.dtype([(f"b{i}",np.uint8) for i in range(n+8)])).reshape(-1)
    view.sort()
    for start in range(1,len(out),8192):
        stop = min(len(out),start+8192)
        assert not np.any(np.all(out[start:stop,:n] == out[start-1:stop-1,:n],axis=1))
    digest = hashlib.sha256(memoryview(out).cast("B")).hexdigest()
    return digest

