"""m_xgdes (abcopy) vs the Fortran original: whole-state comparison."""
from ctypes import byref, c_int

import numpy as np
import pytest

from conftest import globals_diff
import test_xfoil as XF


@pytest.mark.parametrize('seed', range(10))
def test_abcopy_matches_fortran(clib, flib, seed):
    rng = np.random.default_rng(200 + seed)
    XF.init_both(clib, flib, rng)
    x, y = XF.buffer_shape(rng, bool(seed % 2), seed % 4)       # doubled points: LE and 0-3 others
    if seed % 5 == 4:                                           # tripled point (leaves a doubled one)
        j = len(x) // 3
        x, y = np.insert(x, j, [x[j], x[j]]), np.insert(y, j, [y[j], y[j]])
    XF.set_both(clib, flib, 'xb', x); XF.set_both(clib, flib, 'yb', y)
    XF.set_both(clib, flib, 'nb', len(x), np.int32)
    XF.set_both(clib, flib, 'n', int(rng.choice([0, len(x)])), np.int32)   # lblini reset or not
    XF.set_both(clib, flib, 'lblini', 1, np.int32)
    XF.set_both(clib, flib, 'lbflap', int(seed % 3 == 0), np.int32)
    XF.set_both(clib, flib, 'xbf', 0.6)
    assert globals_diff(clib, flib) == []
    clib.abcopy(c_int(1)); flib.__m_xgdes_MOD_abcopy(byref(c_int(1)))
    assert globals_diff(clib, flib) == []


@pytest.mark.parametrize('nb', [0, 1, 640 - 4])
def test_abcopy_rejects_buffer(clib, flib, nb):
    XF.init_both(clib, flib, np.random.default_rng(nb))
    XF.set_both(clib, flib, 'nb', nb, np.int32)
    clib.abcopy(c_int(1)); flib.__m_xgdes_MOD_abcopy(byref(c_int(1)))
    assert globals_diff(clib, flib) == []
