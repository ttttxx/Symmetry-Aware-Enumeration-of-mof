"""Quota prefixes, exact stabilizer filters and certificate scope against brute force."""
import contextlib
import io
import itertools
from collections import Counter

import pytest

from enumerator_rho_quota import RhoQuotaMOFEnumerator
from enumerator_rho_native import HybridRhoQuotaMOFEnumerator
import enumerator_native as native
from tests.test_subgroup_task_dedup import dihedral


def run(cls, perms, quota=3, labels=None, **kw):
    labels = labels or ['Z', 'A']
    e = cls(perms,len(perms[0]),verbose=False,rho_quota_per_level=quota,
            max_total_results=kw.pop('max_total_results',None),**kw)
    with contextlib.redirect_stdout(io.StringIO()):
        records=e.enumerate(labels)
    return e,records


@pytest.mark.parametrize('backend',['python','cpp'])
@pytest.mark.parametrize('prefix',[False,True])
@pytest.mark.parametrize('quota',[1,3,1000])
@pytest.mark.parametrize('perms',[dihedral(4),dihedral(5),list(itertools.permutations(range(3)))])
def test_exact_layer_prefix_matches_original_and_oracle(backend,prefix,quota,perms):
    if backend=='cpp' and not native.native_available(): pytest.skip('build native core')
    old,expected=run(RhoQuotaMOFEnumerator,perms,quota)
    new,actual=run(HybridRhoQuotaMOFEnumerator,perms,quota,backend=backend,prefix_trie=prefix)
    assert actual==expected
    oracle={}
    for c in itertools.product(['A','Z'],repeat=len(perms[0])):
        images={tuple(c[p.index(i)] for i in range(len(c))) for p in perms}
        if c==min(images): oracle[c]=(len(perms)//len(images),len(images))
    histogram=Counter()
    for r in actual:
        key=tuple(r['coloring']);assert oracle[key]==(r['stabilizer_size'],r['orbit_size'])
        histogram[r['stabilizer_size']]+=1
        assert new.classes[new._trie_lookup(key)] is r
    assert all(count<=quota for count in histogram.values())
    if quota==1000:
        assert len(actual)==len(oracle)
        assert sum(r['orbit_size'] for r in actual)==2**len(perms[0])
    assert [(r['stabilizer_order'],r['returned'],r['status']) for r in new.rho_layer_report]==[
        (r['stabilizer_order'],r['returned'],r['status']) for r in old.rho_layer_report]


@pytest.mark.parametrize('backend',['python','cpp'])
@pytest.mark.parametrize('cap',[0,1,3,5,8])
def test_global_cap_and_reuse_reset(backend,cap):
    if backend=='cpp' and not native.native_available():pytest.skip('build native core')
    _,expected=run(RhoQuotaMOFEnumerator,dihedral(4),max_total_results=cap)
    new,actual=run(HybridRhoQuotaMOFEnumerator,dihedral(4),backend=backend,
                   prefix_trie=True,max_total_results=cap)
    assert actual==expected
    with contextlib.redirect_stdout(io.StringIO()):
        assert new.enumerate(['Z','A'])==expected
    assert len(new.classes)==len(expected)
    if cap==0:assert all(r['status']=='not_started_global_output_cap' or r['status']=='global_output_cap' for r in new.rho_layer_report)


@pytest.mark.parametrize('backend',['python','cpp'])
def test_expired_tasks_are_not_layer_coverage(backend):
    if backend=='cpp' and not native.native_available():pytest.skip('build native core')
    e,r=run(HybridRhoQuotaMOFEnumerator,dihedral(4),backend=backend,
             prefix_trie=True,task_time_limit=1e-12)
    assert r==[]
    assert all(x['status']=='time_limited' for x in e.rho_layer_report)
    assert not e._covered_partitions
    assert e.task_dedup_stats['completed_tasks']==0


def test_legacy_quota_task_deadline_is_passed_to_deadline_argument():
    class Capture(RhoQuotaMOFEnumerator):
        def _dfs_lattice_pruning(self,*args,**kw):
            assert len(args)==5
            assert kw['deadline'] is not None
    run(Capture,dihedral(4),task_time_limit=1)


def test_factory_defaults_to_optimized_auto_and_selects_quota_backend():
    from main import get_enumerator_implementation
    from enumerator import MOFEnumerator
    from enumerator_native import HybridMOFEnumerator
    assert get_enumerator_implementation()==(HybridMOFEnumerator,{'backend':'auto'})
    assert get_enumerator_implementation(rho_quota_per_level=2)==(HybridRhoQuotaMOFEnumerator,{'backend':'auto','rho_quota_per_level':2})
    assert get_enumerator_implementation('legacy')==(MOFEnumerator,{})
    assert get_enumerator_implementation('legacy',rho_quota_per_level=2)==(RhoQuotaMOFEnumerator,{'rho_quota_per_level':2})
    for backend in ['python','cpp','auto']:
        cls,opts=get_enumerator_implementation(backend,True,0,'trie',2)
        assert cls is HybridRhoQuotaMOFEnumerator
        assert opts['backend']==backend and opts['rho_quota_per_level']==2


def test_auto_falls_back_to_optimized_python(monkeypatch):
    def absent():raise OSError('missing test library')
    monkeypatch.setattr(native,'_load_core',absent)
    _,expected=run(RhoQuotaMOFEnumerator,dihedral(4))
    e,r=run(HybridRhoQuotaMOFEnumerator,dihedral(4),backend='auto',prefix_trie=True)
    assert r==expected and e.active_backend=='python'


@pytest.mark.parametrize('store',['hash','trie'])
def test_native_quota_store_and_exact_layer_scope(store):
    if not native.native_available():pytest.skip('build native core')
    _,expected=run(RhoQuotaMOFEnumerator,dihedral(5),quota=1000,labels=['C','A','B'])
    e,r=run(HybridRhoQuotaMOFEnumerator,dihedral(5),quota=1000,labels=['C','A','B'],
             backend='cpp',native_store=store,prefix_trie=True,prefix_trie_max_nodes=0)
    assert r==expected
    assert len({tuple(x['coloring']) for x in r})==len(r)
    assert any(x['stabilizer_size']==1 for x in r)
    assert e._coverage_order==1
