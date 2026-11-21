"""api.f90 entry points (C src/xfoil_api.c vs the Fortran bind(c) originals): call sequences
through the public API, whole-state comparison after every call, getters compared by value."""
from ctypes import byref, c_bool, c_double, c_int, c_void_p

import numpy as np
import pytest

from conftest import global_array, globals_diff, sync_globals
import test_xfoil as XF

GETTERS = {'get_n_coords': c_int, 'get_reynolds': c_double, 'get_mach': c_double,
           'get_n_crit': c_double, 'get_max_iter': c_int, 'get_print': c_bool}


def call(clib, flib, name, *args):
    """Call `name` in both libraries with fresh by-reference copies of the scalar args;
    numpy arrays are copied per library. Returns the two results."""
    out = []
    for lib in (clib, flib):
        conv = []
        for a in args:
            if isinstance(a, np.ndarray):
                conv.append(np.array(a, dtype=np.float64).ctypes.data_as(c_void_p))
            elif isinstance(a, (int, np.integer)):
                conv.append(byref(c_int(int(a))))
            else:
                conv.append(byref(c_double(float(a))))
        f = getattr(lib, name)
        f.restype = GETTERS.get(name)
        out.append(f(*conv))
    assert globals_diff(clib, flib) == [], name
    return out


def check_getters(clib, flib):
    for name in GETTERS:
        c, f = call(clib, flib, name)
        assert c == f, name
    vals = []
    for lib in (clib, flib):
        t, b = c_double(), c_double()
        lib.get_xtr(byref(t), byref(b))
        vals.append((t.value, b.value))
    assert vals[0] == vals[1]
    n = clib.get_n_coords()
    xs = []
    for lib in (clib, flib):
        x, y = np.zeros(n), np.zeros(n)
        lib.get_airfoil(x.ctypes.data_as(c_void_p), y.ctypes.data_as(c_void_p), byref(c_int(n)))
        xs.append(np.r_[x, y].tobytes())
    assert xs[0] == xs[1]


def fresh(clib, flib):
    sync_globals(clib, flib)
    clib.init(); flib.init()
    assert globals_diff(clib, flib) == []


@pytest.mark.parametrize('seed', range(12))
def test_api_sequence_matches_fortran(clib, flib, seed):
    rng = np.random.default_rng(400 + seed)
    fresh(clib, flib)
    x, y = XF.buffer_shape(rng, bool(seed % 2), seed % 3)
    if seed % 4 == 1:                                  # clockwise input -> reversed
        x, y = x[::-1].copy(), y[::-1].copy()
    if seed % 6 == 2:
        XF.set_both(clib, flib, 'lnorm', 1, np.int32)
    call(clib, flib, 'set_airfoil', x, y, len(x))
    check_getters(clib, flib)
    if seed % 3 == 0:
        call(clib, flib, 'set_naca', int(rng.choice([12, 2412, 4415, 23012, 0])))
    call(clib, flib, 'repanel', int(rng.choice([60, 160, 300, 700])), rng.uniform(0.5, 2), rng.uniform(0.05, 0.5),
         rng.uniform(0.1, 0.5), 1.0, 1.0, *sorted(rng.uniform(0.1, 0.9, 2)))
    call(clib, flib, 'set_reynolds', float(rng.choice([0.0, 1e6, 3e5])))
    for m in (0.0, float(rng.uniform(0.1, 0.6)), 0.3, 0.3):
        call(clib, flib, 'set_mach', m)
    call(clib, flib, 'set_xtr', rng.uniform(0, 1), rng.uniform(0, 1))
    call(clib, flib, 'set_n_crit', rng.uniform(4, 12))
    call(clib, flib, 'set_max_iter', int(rng.integers(1, 200)))
    call(clib, flib, 'reset_bls')
    check_getters(clib, flib)


def test_get_print_and_set_print(clib, flib):
    for v in (True, False):
        for lib in (clib, flib):
            lib.set_print(byref(c_bool(v)))
            lib.get_print.restype = c_bool
            assert lib.get_print() is v
    assert globals_diff(clib, flib) == []
