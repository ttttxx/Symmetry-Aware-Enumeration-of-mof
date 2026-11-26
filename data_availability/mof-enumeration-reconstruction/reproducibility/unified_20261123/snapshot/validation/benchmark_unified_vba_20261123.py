"""Current C++ bounded-output benchmarks using the common lexical label pools."""
import argparse
import collections
import contextlib
import gc
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def worker(c4, policy):
    import sympy.combinatorics
    from enumerator_native import HybridMOFEnumerator
    from enumerator_rho_native import HybridRhoQuotaMOFEnumerator
    from experiment.casm_scale_common import pools_for, orbit_count
    from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
    from symmetry import SymmetryEngine
    pool = pools_for(1, c4)[4]
    cap = {2: 41984, 3: 60549}[c4]
    with contextlib.redirect_stdout(sys.stderr):
        graph = enrich_graph_data_with_topology(load_graph(str(ROOT/'data/vba/vba_graph.json')))
        degrees = [node['degree'] for node in graph['nodes']]
        perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
        cls = HybridRhoQuotaMOFEnumerator if policy == 'quota' else HybridMOFEnumerator
        extra = dict(rho_quota_per_level=20000) if policy == 'quota' else dict(identity_orbit_limit=float('inf'))
        instance = cls(perms, len(degrees), verbose=False, backend='cpp',
            native_store='trie', prefix_trie=True, user_selected_ligands={4: pool},
            topology_requirements=calculate_topology_requirements(graph), node_degrees=degrees,
            max_total_results=100000 if policy == 'quota' else cap, task_time_limit=60, **extra)
        gc.collect()
        start = time.perf_counter()
        records = instance.enumerate(sorted(pool))
        elapsed = time.perf_counter()-start
    assert instance.active_backend == 'cpp' and instance.active_store == 'cpp_trie'
    assert len(records) == cap
    distribution = {}
    for order in [16, 8, 4, 2, 1]:
        subset = [row for row in records if row['stabilizer_size'] == order]
        keys = {tuple(row['coloring']) for row in subset}
        assert len(keys) == len(subset)
        assert all(row['orbit_size'] * order == len(perms) for row in subset)
        vectors = {tuple(collections.Counter(row['coloring'])[label] for label in pool) for row in subset}
        distribution[order] = dict(returned=len(subset), composition_vectors=len(vectors))
    expected = ({2: [4,90,1890,20000,20000], 3: [9,540,20000,20000,20000]} if policy == 'quota'
                else {2: [4,90,1890,40000,0], 3: [9,540,55440,4560,0]})[c4]
    assert [distribution[s]['returned'] for s in [16,8,4,2,1]] == expected
    result = dict(c4=c4, policy=policy, pool=pool, python=sys.version, native_abi=4,
        group_size=len(perms), naive=c4**len(degrees), returned=len(records), time_s=elapsed,
        burnside_total=orbit_count(perms,degrees,{4:pool}),
        distribution=distribution, backend=instance.active_backend, store=instance.active_store,
        ordered_output_sha256=hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest(),
        peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1048576 if sys.platform=='darwin' else 1024),
        native_kernel_s=sum(row['kernel_s'] for row in instance.native_task_reports),
        prefix_stats=instance.prefix_trie_stats, task_stats=instance.task_dedup_stats,
        output_cap=100000 if policy=='quota' else cap, task_time_limit=60)
    if policy == 'quota':
        result['layers'] = instance.rho_layer_report
        assert all(row['status'] in ('quota_reached','layer_exhausted') for row in result['layers'])
    else:
        assert not any(row['timed_out'] for row in instance.native_task_reports)
    return result


def summarize(records):
    summary = []
    for c4 in [2,3]:
        row = dict(c4=c4)
        for policy in ['prefix','quota']:
            runs = [r for r in records if r['c4']==c4 and r['policy']==policy]
            assert len(runs)==3 and len({r['ordered_output_sha256'] for r in runs})==1
            assert all(r['distribution']==runs[0]['distribution'] for r in runs)
            item = dict(pool=runs[0]['pool'],returned=runs[0]['returned'],distribution=runs[0]['distribution'],
                median_s=statistics.median(r['time_s'] for r in runs),
                min_s=min(r['time_s'] for r in runs),max_s=max(r['time_s'] for r in runs),
                ordered_output_sha256=runs[0]['ordered_output_sha256'])
            if policy=='quota':
                item['layers'] = [dict(stabilizer_order=s,
                    median_s=statistics.median(next(x['elapsed_s'] for x in r['layers'] if x['stabilizer_order']==s) for r in runs),
                    status=next(x['status'] for x in runs[0]['layers'] if x['stabilizer_order']==s)) for s in [16,8,4,2,1]]
            row[policy]=item
        assert row['prefix']['returned']==row['quota']['returned']
        summary.append(row)
    return summary


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker',type=int,choices=[2,3])
    parser.add_argument('--policy',choices=['prefix','quota'])
    parser.add_argument('--output',type=Path,default=ROOT/'validation/MATCH_UNIFIED_VBA_2026-11-23.json')
    args=parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker,args.policy)))
        sys.exit()
    checkpoint=args.output.with_suffix('.jsonl')
    assert not args.output.exists() and not checkpoint.exists(), 'Use a fresh output path'
    rows=[]
    for c4 in [2,3]:
        for rep in range(3):
            for policy in (['prefix','quota'] if rep%2==0 else ['quota','prefix']):
                done=subprocess.run([sys.executable,__file__,'--worker',str(c4),'--policy',policy],
                    capture_output=True,text=True,timeout=180)
                if done.returncode:
                    raise RuntimeError(done.stderr)
                row=json.loads(done.stdout)
                row['rep']=rep
                rows.append(row)
                with checkpoint.open('a') as stream:
                    stream.write(json.dumps(row)+'\n')
                print(c4,policy,rep+1,f"{row['time_s']:.4f}s",flush=True)
    args.output.write_text(json.dumps(dict(records=rows,summary=summarize(rows),platform=platform.platform(),
        native_abi=4,timing_scope='enumerate; imports/group construction and validation excluded; task setup and outputs included'),indent=2)+'\n')
