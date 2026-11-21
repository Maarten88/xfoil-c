# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""m_spline.f90 -> src/xspline.c: C results must be bit-identical to the Fortran routines."""
import os
import subprocess
import sys
import textwrap
from ctypes import byref, c_double, c_int, c_void_p

import numpy as np
import pytest

from conftest import FORTRAN_LIB

N_CASES = 150
FUNCS = ('seval', 'deval', 'd2val', 'curv')


def ptr(arr):
    return arr.ctypes.data_as(c_void_p)   # typed pointer: a plain int is passed as 32-bit C int


@pytest.fixture(scope='module')
def libs(clib, flib):
    for name in FUNCS:
        getattr(clib, name).restype = c_double
        getattr(flib, '__m_spline_MOD_' + name).restype = c_double
    return clib, flib


def F(flib, name):
    return getattr(flib, '__m_spline_MOD_' + name)


def cases(n_cases=N_CASES):
    for seed in range(n_cases):
        yield np.random.default_rng(1000 + seed)


def random_curve(rng, n=None):
    """Strictly increasing s with uneven spacing, and smooth-ish x(s), y(s)."""
    n = int(rng.integers(2, 300)) if n is None else n
    s = np.cumsum(rng.uniform(0.01, 1.0, size=n))
    s -= s[0]
    x = np.sin(s * rng.uniform(0.2, 2.0)) + 0.1 * rng.normal(size=n)
    y = np.cos(s * rng.uniform(0.2, 2.0)) + 0.1 * rng.normal(size=n)
    return s, x, y


def segmented_curve(rng):
    """s with segment joints (equal successive values), each segment >= 2 points."""
    lengths = rng.integers(2, 40, size=int(rng.integers(1, 6)))
    s_parts, s0 = [], 0.0
    for m in lengths:
        seg = s0 + np.concatenate([[0.0], np.cumsum(rng.uniform(0.01, 1.0, size=m - 1))])
        s_parts.append(seg)
        s0 = seg[-1]                    # next segment starts at the same s -> joint
    s = np.concatenate(s_parts)
    x = np.sin(s) + 0.1 * rng.normal(size=len(s))
    return s, x


def run_spline_c(clib, name, x, s, *extra):
    xs = np.zeros_like(x)
    getattr(clib, name)(ptr(x), ptr(xs), ptr(s), c_int(len(x)), *[c_double(v) for v in extra])
    return xs


def run_spline_f(flib, name, x, s, *extra):
    xs = np.zeros_like(x)
    F(flib, name)(ptr(x), ptr(xs), ptr(s), byref(c_int(len(x))), *[byref(c_double(v)) for v in extra])
    return xs


# ------------------------------------------------------------------ coefficient routines

@pytest.mark.parametrize('xs1', [999.0, -999.0, 0.37])
@pytest.mark.parametrize('xs2', [999.0, -999.0, -1.25])
def test_splind_matches_fortran(libs, xs1, xs2):
    clib, flib = libs
    for rng in cases(40):
        s, x, _ = random_curve(rng)
        c = run_spline_c(clib, 'splind', x, s, xs1, xs2)
        assert np.array_equal(c, run_spline_f(flib, 'splind', x, s, xs1, xs2))


def test_splind_two_points_zero_third_derivative(libs):
    clib, flib = libs
    s, x = np.array([0.0, 0.7]), np.array([1.0, -0.4])
    c = run_spline_c(clib, 'splind', x, s, -999.0, -999.0)
    assert np.array_equal(c, run_spline_f(flib, 'splind', x, s, -999.0, -999.0))


def test_trisol_matches_fortran(libs):
    clib, flib = libs
    for rng in cases():
        kk = int(rng.integers(1, 200))
        a = rng.uniform(2.0, 4.0, kk)
        b, c, d = rng.normal(size=kk), rng.normal(size=kk), rng.normal(size=kk)
        ca, cc, cd = a.copy(), c.copy(), d.copy()
        fa, fc, fd = a.copy(), c.copy(), d.copy()
        clib.trisol(ptr(ca), ptr(b), ptr(cc), ptr(cd), c_int(kk))
        F(flib, 'trisol')(ptr(fa), ptr(b), ptr(fc), ptr(fd), byref(c_int(kk)))
        assert np.array_equal(cd, fd) and np.array_equal(ca, fa) and np.array_equal(cc, fc)


# ------------------------------------------------------------------ evaluation routines

def eval_points(rng, s):
    """Random points inside and slightly outside the range, plus every knot."""
    span = s[-1] - s[0]
    return np.concatenate([rng.uniform(s[0] - 0.1 * span, s[-1] + 0.1 * span, 30), s])


def test_seval_deval_d2val_match_fortran(libs):
    clib, flib = libs
    for rng in cases(60):
        s, x, _ = random_curve(rng)
        xs = run_spline_c(clib, 'segspl', x, s)
        n = len(s)
        for ss in eval_points(rng, s):
            for name in ('seval', 'deval', 'd2val'):
                c = getattr(clib, name)(c_double(ss), ptr(x), ptr(xs), ptr(s), c_int(n))
                f = F(flib, name)(byref(c_double(ss)), ptr(x), ptr(xs), ptr(s), byref(c_int(n)))
                assert c == f, (name, ss)


def test_curv_matches_fortran(libs):
    clib, flib = libs
    for rng in cases(60):
        s, x, y = random_curve(rng)
        xs = run_spline_c(clib, 'segspl', x, s)
        ys = run_spline_c(clib, 'segspl', y, s)
        n = len(s)
        for ss in eval_points(rng, s):
            for name in ('curv',):
                c = getattr(clib, name)(c_double(ss), ptr(x), ptr(xs), ptr(y), ptr(ys), ptr(s), c_int(n))
                f = F(flib, name)(byref(c_double(ss)), ptr(x), ptr(xs), ptr(y), ptr(ys), ptr(s),
                                  byref(c_int(n)))
                assert c == f, (name, ss)


def test_sinvrt_matches_fortran(libs):
    clib, flib = libs
    converged = 0
    for rng in cases():
        n = int(rng.integers(3, 100))
        s = np.cumsum(rng.uniform(0.05, 1.0, size=n)) - 0.05
        x = s + 0.2 * np.sin(s)          # monotonic, so the inverse exists
        xs = run_spline_c(clib, 'segspl', x, s)
        xi = float(rng.uniform(x[0], x[-1]))
        guess = float(rng.uniform(s[0], s[-1]))
        c_si, f_si = c_double(guess), c_double(guess)
        clib.sinvrt(byref(c_si), c_double(xi), ptr(x), ptr(xs), ptr(s), c_int(n))
        F(flib, 'sinvrt')(byref(f_si), byref(c_double(xi)), ptr(x), ptr(xs), ptr(s), byref(c_int(n)))
        assert c_si.value == f_si.value
        converged += c_si.value != guess
    assert converged > N_CASES * 0.8


def test_scalc_matches_fortran(libs):
    clib, flib = libs
    for rng in cases():
        _, x, y = random_curve(rng)
        cs, fs = np.zeros_like(x), np.zeros_like(x)
        clib.scalc(ptr(x), ptr(y), ptr(cs), c_int(len(x)))
        F(flib, 'scalc')(ptr(x), ptr(y), ptr(fs), byref(c_int(len(x))))
        assert np.array_equal(cs, fs)


def test_segspl_matches_fortran(libs):
    clib, flib = libs
    joints = 0
    for rng in cases():
        s, x = segmented_curve(rng)
        joints += int(np.sum(np.diff(s) == 0))
        assert np.array_equal(run_spline_c(clib, 'segspl', x, s), run_spline_f(flib, 'segspl', x, s))
    assert joints > N_CASES          # inputs really contain segment joints


# ------------------------------------------------------------------ fatal errors

STOP_SCRIPT = textwrap.dedent('''
    import sys
    from ctypes import CDLL, byref, c_int, c_void_p
    import numpy as np
    lib = CDLL(sys.argv[1])
    s = np.array([0.0, 0.0, 1.0, 2.0]); x = np.zeros(4); xs = np.zeros(4)
    if sys.argv[2] == 'c':
        lib.segspl(x.ctypes.data_as(c_void_p), xs.ctypes.data_as(c_void_p), s.ctypes.data_as(c_void_p), c_int(4))
    else:
        lib.__m_spline_MOD_segspl(x.ctypes.data_as(c_void_p), xs.ctypes.data_as(c_void_p), s.ctypes.data_as(c_void_p), byref(c_int(4)))
    print('NOT STOPPED')
''')


def test_stop_matches_fortran(clib, flib):
    from xfoil_c import _lib
    results = {}
    for tag, path in (('c', _lib.lib_path), ('f', FORTRAN_LIB)):
        p = subprocess.run([sys.executable, '-c', STOP_SCRIPT, path, tag], capture_output=True, text=True)
        results[tag] = (p.returncode, p.stdout.strip(), p.stderr.strip())
    assert results['c'] == results['f']
    assert results['c'] == (0, '', 'STOP SEGSPL:  First input point duplicated')
