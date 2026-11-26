"""Independent checks for the scale benchmark's coverage and normalization."""
import itertools

import pytest

from experiment.casm_scale_common import orbit_count, record_digest


def test_cycle_count_matches_independent_brute_force():
    perms=[(0,1,2),(1,2,0),(2,0,1)]
    classes={min(tuple(c[p[i]] for i in range(3)) for p in perms)
             for c in itertools.product("AB",repeat=3)}
    assert orbit_count(perms,[3]*3,{3:["A","B"]}) == len(classes) == 4


def test_record_digest_normalizes_canonical_convention_and_order():
    perms=[(0,1,2),(1,2,0),(2,0,1)]
    a=[dict(id=0,coloring=list("AAA"),stabilizer_size=3,orbit_size=1),
       dict(id=1,coloring=list("AAB"),stabilizer_size=1,orbit_size=3)]
    b=[dict(id=0,coloring=list("BAA"),stabilizer_size=1,orbit_size=3),
       dict(id=1,coloring=list("AAA"),stabilizer_size=3,orbit_size=1)]
    assert record_digest(a,["A","B"],perms) == record_digest(b,["A","B"],perms)
    b[0]["stabilizer_size"]=3
    assert record_digest(a,["A","B"],perms) != record_digest(b,["A","B"],perms)


def test_record_digest_rejects_two_representatives_of_same_orbit():
    rows=[dict(coloring=list(c),stabilizer_size=1,orbit_size=3) for c in ["AAB","BAA"]]
    with pytest.raises(AssertionError):
        record_digest(rows,["A","B"],[(0,1,2),(1,2,0),(2,0,1)])
