"""Fresh-process exact-rho quota validation and timing on the paper's vba pools."""
import argparse,collections,contextlib,gc,hashlib,json,os
from pathlib import Path
import platform,resource,statistics,subprocess,sys,time
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
SOURCES=['enumerator.py','enumerator_native.py','enumerator_rho_native.py','enumerator_rho_quota.py',
'enumerator_casm_inspired.py','enumerator_task_dedup.py','enumerator_storage.py','search_prefix_trie.py',
'native/enumeration_core.cpp','native/libmof_enumeration.dylib','main.py','symmetry.py',
'validation/benchmark_rho_native_20261123.py','data/vba/vba_graph.json']

def worker(c4,backend):
    import sympy.combinatorics
    from enumerator_rho_quota import RhoQuotaMOFEnumerator
    from enumerator_rho_native import HybridRhoQuotaMOFEnumerator
    from main import load_graph,enrich_graph_data_with_topology,calculate_topology_requirements
    from symmetry import SymmetryEngine
    pool=['4c_Cu_1_Ch','4c_B_Ch']+(['4c_Cd_1_Ch'] if c4==3 else [])
    with contextlib.redirect_stdout(sys.stderr):
        graph=enrich_graph_data_with_topology(load_graph(str(ROOT/'data/vba/vba_graph.json')))
        perms=SymmetryEngine(graph,verbose=False).compute_automorphism_group()
        cls=RhoQuotaMOFEnumerator if backend=='legacy' else HybridRhoQuotaMOFEnumerator
        extra={} if backend=='legacy' else dict(backend=backend,native_store='trie',prefix_trie=True)
        e=cls(perms,len(graph['nodes']),verbose=False,user_selected_ligands={4:pool},
              topology_requirements=calculate_topology_requirements(graph),
              node_degrees=[n['degree'] for n in graph['nodes']],rho_quota_per_level=20000,
              max_total_results=100000,task_time_limit=60,**extra)
        gc.collect();start=time.perf_counter();records=e.enumerate(sorted(pool));elapsed=time.perf_counter()-start
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1048576 if sys.platform=='darwin' else 1024)
    distribution={}
    for order in [16,8,4,2,1]:
        rr=[r for r in records if r['stabilizer_size']==order]
        vectors={tuple(collections.Counter(r['coloring'])[label] for label in pool) for r in rr}
        distribution[order]=dict(returned=len(rr),composition_vectors=len(vectors))
    assert [distribution[s]['returned'] for s in [16,8,4,2,1]]==({2:[4,90,1890,20000,20000],3:[9,540,20000,20000,20000]}[c4])
    assert all(r['status'] in ('quota_reached','layer_exhausted') for r in e.rho_layer_report)
    digest=hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest()
    result=dict(c4=c4,backend=backend,python=sys.version,pool=pool,group_size=len(perms),
        naive=c4**len(graph['nodes']),returned=len(records),time_s=elapsed,peak_rss_mib=rss,
        ordered_output_sha256=digest,distribution=distribution,layers=e.rho_layer_report)
    if backend!='legacy':
        assert e.active_backend==backend
        result.update(task_stats=e.task_dedup_stats,prefix_stats=e.prefix_trie_stats,
            native_store=e.active_store,index_stats=e.result_index_stats,
            native_kernel_s=sum(r['kernel_s'] for r in e.native_task_reports))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker',type=int,choices=[2,3]);parser.add_argument('--backend',choices=['legacy','python','cpp'])
    parser.add_argument('--reps',type=int,default=3);parser.add_argument('--output',type=Path,default=ROOT/'validation/RHO_NATIVE_VBA_2026-11-23.json')
    args=parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker,args.backend)));sys.exit()
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}
    rows=[]
    for c4 in [2,3]:
        for rep in range(args.reps):
            sides=['legacy','python','cpp'];offset=rep%3
            for side in sides[offset:]+sides[:offset]:
                r=subprocess.run([sys.executable,__file__,'--worker',str(c4),'--backend',side],capture_output=True,text=True,timeout=180)
                if r.returncode:raise RuntimeError(r.stderr)
                row=json.loads(r.stdout);row['rep']=rep;rows.append(row)
                with args.output.with_suffix('.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
                print(c4,side,rep,f"{row['time_s']:.4f}s",row['returned'],flush=True)
    summary=[]
    for c4 in [2,3]:
        rr=[r for r in rows if r['c4']==c4]
        assert len({r['ordered_output_sha256'] for r in rr})==1
        s=dict(c4=c4,returned=rr[0]['returned'],distribution=rr[0]['distribution'],all_ordered_outputs_match=True)
        for side in ['legacy','python','cpp']:
            x=[r for r in rr if r['backend']==side]
            s[side]=dict(median_s=statistics.median(r['time_s'] for r in x),
                         median_rss_mib=statistics.median(r['peak_rss_mib'] for r in x))
        summary.append(s)
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in hashes.items())
    args.output.write_text(json.dumps(dict(records=rows,summary=summary,source_sha256=hashes,
        platform=platform.platform(),native_abi=4,timing_scope='enumerate; imports/group construction and verification excluded; task preparation/output included'),indent=2)+'\n')
