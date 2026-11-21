# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""m_xutils.f90 -> src/xutils.c: C results must be bit-identical to the Fortran routines."""
import math
import subprocess
import sys
import textwrap
from ctypes import byref, c_double, c_int, c_void_p

import numpy as np
import pytest

from conftest import FORTRAN_LIB

N_CASES = 300


def ptr(arr):
    return arr.ctypes.data_as(c_void_p)   # typed pointer: a plain int is passed as 32-bit C int


@pytest.fixture(scope='module')
def libs(clib, flib):
    clib.atanc.restype = c_double
    flib.__m_xutils_MOD_atanc.restype = c_double
    return clib, flib


def run_setexp(clib, flib, ds1, smax, nn):
    cs, fs = np.zeros(nn), np.zeros(nn)
    clib.setexp(ptr(cs), c_double(ds1), c_double(smax), c_int(nn))
    flib.__m_xutils_MOD_setexp(ptr(fs), byref(c_double(ds1)), byref(c_double(smax)), byref(c_int(nn)))
    return cs, fs


def test_setexp_matches_fortran(libs):
    clib, flib = libs
    kinds = {'stretched': 0, 'compressed': 0, 'uniform': 0}
    for seed in range(N_CASES):
        rng = np.random.default_rng(seed)
        nn = int(rng.integers(3, 300))
        ds1 = float(rng.uniform(1e-4, 0.1))
        kind = ('stretched', 'compressed', 'uniform')[seed % 3]
        if kind == 'stretched':          # ratio > 1
            smax = (nn - 1) * ds1 * float(rng.uniform(1.01, 50.0))
        elif kind == 'compressed':       # ratio < 1
            smax = (nn - 1) * ds1 * float(rng.uniform(0.3, 0.99))
        else:                            # ratio guess exactly 1: Newton skipped
            smax = (nn - 1) * ds1
        kinds[kind] += 1
        cs, fs = run_setexp(clib, flib, ds1, smax, nn)
        assert cs.tobytes() == fs.tobytes(), (seed, kind, nn)
    assert min(kinds.values()) > 0


def test_setexp_three_points(libs):
    clib, flib = libs                    # nex == 2: closed-form ratio
    for smax in (0.05, 0.3, 1.0, 7.0):
        cs, fs = run_setexp(clib, flib, 0.1, smax, 3)
        assert cs.tobytes() == fs.tobytes()


def test_setexp_nonconvergence_matches_fortran(libs):
    clib, flib = libs                    # extreme stretching: Newton hits 100 iterations
    cs, fs = run_setexp(clib, flib, 1e-12, 1e6, 400)
    assert cs.tobytes() == fs.tobytes()


def atanc_both(clib, flib, y, x, th):
    c = clib.atanc(c_double(y), c_double(x), c_double(th))
    f = flib.__m_xutils_MOD_atanc(byref(c_double(y)), byref(c_double(x)), byref(c_double(th)))
    return c, f


def test_atanc_matches_fortran(libs):
    clib, flib = libs
    rng = np.random.default_rng(7)
    wraps = 0
    for _ in range(20000):
        r = float(rng.uniform(0.01, 10.0))
        a = float(rng.uniform(-math.pi, math.pi))
        y, x = r * math.sin(a), r * math.cos(a)
        th = a + float(rng.uniform(-3.0, 3.0)) + 2 * math.pi * int(rng.integers(-10, 11))
        c, f = atanc_both(clib, flib, y, x, th)
        assert c.hex() == f.hex(), (y, x, th)
        wraps += abs(c - math.atan2(y, x)) > 1.0
    assert wraps > 10000                 # most cases needed a 2*pi correction


@pytest.mark.parametrize('y,x,th', [
    (-1.0, -1.0, 0.75 * math.pi),       # docstring example: 1.25*pi
    (0.0, 1.0, 0.0), (-0.0, 1.0, 0.0),  # dthet = +0.0 / -0.0
    (0.0, -1.0, math.pi), (-0.0, -1.0, math.pi), (0.0, -1.0, -math.pi),
    (1e-300, -1.0, 3.0), (0.0, 0.0, 0.0), (1.0, 0.0, 1e6), (1.0, 0.0, -1e9),
])
def test_atanc_edge_cases(libs, y, x, th):
    c, f = atanc_both(*libs, y, x, th)
    assert c.hex() == f.hex()


def test_atanc_docstring_example(libs):
    c, _ = atanc_both(*libs, -1.0, -1.0, 0.75 * math.pi)
    assert c == pytest.approx(1.25 * math.pi)


STOP_SCRIPT = textwrap.dedent('''
    import sys
    from ctypes import CDLL, byref, c_double, c_int, c_void_p
    import numpy as np
    lib = CDLL(sys.argv[1]); s = np.zeros(2)
    if sys.argv[2] == 'c':
        lib.setexp(s.ctypes.data_as(c_void_p), c_double(0.1), c_double(1.0), c_int(2))
    else:
        lib.__m_xutils_MOD_setexp(s.ctypes.data_as(c_void_p), byref(c_double(0.1)), byref(c_double(1.0)), byref(c_int(2)))
    print('NOT STOPPED')
''')


def test_setexp_stop_matches_fortran(clib, flib):
    from xfoil_c import _lib
    results = {}
    for tag, path in (('c', _lib.lib_path), ('f', FORTRAN_LIB)):
        p = subprocess.run([sys.executable, '-c', STOP_SCRIPT, path, tag], capture_output=True, text=True)
        results[tag] = (p.returncode, p.stdout.strip(), p.stderr.strip())
    assert results['c'] == results['f'] == (0, '', 'STOP SETEXP: Cannot fill array.  N too small.')
