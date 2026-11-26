"""Same-size bounded vba legacy/current comparison; separate from exact-rho quotas."""
import argparse, collections, contextlib, gc, hashlib, json, os
from pathlib import Path
import platform, resource, statistics, subprocess, sys, time
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

def worker(c4,side):
    import sympy.combinatorics
    from enumerator import MOFEnumerator
    from enumerator_native import HybridMOFEnumerator
    from main import load_graph,enrich_graph_data_with_topology,calculate_topology_requirements
    from symmetry import SymmetryEngine
    from experiment.casm_scale_common import orbit_count,record_digest
    pool=['4c_Cu_1_Ch','4c_B_Ch']+(['4c_Cd_1_Ch'] if c4==3 else [])
    cap={2:41984,3:60549}[c4]
    with contextlib.redirect_stdout(sys.stderr):
        graph=enrich_graph_data_with_topology(load_graph(str(ROOT/'data/vba/vba_graph.json')))
        degrees=[node['degree'] for node in graph['nodes']]
        perms=SymmetryEngine(graph,verbose=False).compute_automorphism_group()
        cls=MOFEnumerator if side=='legacy' else HybridMOFEnumerator
        extra={} if side=='legacy' else dict(backend='cpp',native_store='trie',prefix_trie=True)
        instance=cls(perms,len(degrees),verbose=False,user_selected_ligands={4:pool},
            node_degrees=degrees,topology_requirements=calculate_topology_requirements(graph),
            max_total_results=cap,task_time_limit=None,identity_orbit_limit=float('inf'),**extra)
        gc.collect()
        start=time.perf_counter();records=instance.enumerate(sorted(pool));elapsed=time.perf_counter()-start
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1048576 if sys.platform=='darwin' else 1024)
    assert len(records)==cap
    result=dict(c4=c4,side=side,time_s=elapsed,peak_rss_mib=peak,returned=len(records),
        group_size=len(perms),naive=c4**len(degrees),burnside_total=orbit_count(perms,degrees,{4:pool}),
        status='output_cap',timeout_disabled=True,python=sys.version,
        ordered_digest=hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest(),
        record_digest=record_digest(records,sorted(pool),perms),
        stabilizer_histogram=dict(collections.Counter(r['stabilizer_size'] for r in records)),
        represented_assignments=sum(r['orbit_size'] for r in records))
    if side=='ours':
        assert instance.active_backend=='cpp' and instance.active_store=='cpp_trie'
        result.update(kernel_s=sum(r['kernel_s'] for r in instance.native_task_reports),
            prefix_stats=instance.prefix_trie_stats,task_stats=instance.task_dedup_stats)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker',type=int,choices=[2,3]);parser.add_argument('--side',choices=['legacy','ours'])
    args=parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker,args.side)));sys.exit()
    runs=[]
    for c4 in [2,3]:
        for rep in range(3):
            for side in (['legacy','ours'] if rep%2==0 else ['ours','legacy']):
                r=subprocess.run([sys.executable,__file__,'--worker',str(c4),'--side',side],capture_output=True,text=True,timeout=180)
                if r.returncode:raise RuntimeError(r.stderr)
                row=json.loads(r.stdout);row['rep']=rep;runs.append(row)
                print(c4,side,rep,row['time_s'],flush=True)
    summaries=[]
    for c4 in [2,3]:
        rows=[r for r in runs if r['c4']==c4]
        assert len({r['ordered_digest'] for r in rows})==1
        assert len({r['record_digest'] for r in rows})==1
        s=dict(c4=c4,all_ordered_outputs_match=True,returned=rows[0]['returned'],
            naive=rows[0]['naive'],burnside_total=rows[0]['burnside_total'],
            stabilizer_histogram=rows[0]['stabilizer_histogram'])
        for side in ['legacy','ours']:
            rr=[r for r in rows if r['side']==side]
            s[side+'_median_s']=statistics.median(r['time_s'] for r in rr)
            s[side+'_median_rss_mib']=statistics.median(r['peak_rss_mib'] for r in rr)
        s['legacy_speedup']=s['legacy_median_s']/s['ours_median_s'];summaries.append(s)
    Path(ROOT/'validation/MATCH_VBA_REFRESH_2026-11-16.json').write_text(json.dumps(dict(records=runs,summary=summaries,platform=platform.platform()),indent=2)+'\n')
