import sys, os, time
sys.path.insert(0, '/Users/jing/PyCharmMiscProject/mof')
from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator

graph_data = load_graph('/Users/jing/PyCharmMiscProject/mof/data/fjh/fjh_graph.json')
graph_data = enrich_graph_data_with_topology(graph_data)
topo_req = calculate_topology_requirements(graph_data)
node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
ligands = {3: ['3c_B_Ch', '3c_triazine_Ch'], 4: ['4c_Cu_1_Ch', '4c_B_Ch', '4c_Cd_1_Ch']}
all_ligs = sorted(set(l for ls in ligands.values() for l in ls))
perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
en = MOFEnumerator(perms, 28, verbose=True,
                   topology_requirements=topo_req,
                   user_selected_ligands=ligands, node_degrees=node_degrees)
t0 = time.time()
results = en.enumerate(all_ligs)
elapsed = time.time() - t0
print(f'fjh large: {len(results)} structures in {elapsed:.2f}s')
assert elapsed < 60, f'Should complete in <60s, took {elapsed:.1f}s'
print('PASS: completed in time')
