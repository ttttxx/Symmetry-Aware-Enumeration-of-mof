import os
import sys
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from enumerator import IDENTITY_ORBIT_LIMIT, IDENTITY_NEW_LIMIT

def test_constants_exist_and_positive():
    assert IDENTITY_ORBIT_LIMIT > 0
    assert IDENTITY_NEW_LIMIT > 0
    assert IDENTITY_ORBIT_LIMIT > IDENTITY_NEW_LIMIT
