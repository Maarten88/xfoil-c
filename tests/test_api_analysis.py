"""api.f90 analysis entry points (alfa, cl, aseq, cseq, get_n_cp, get_cp): C vs Fortran end to end,
on the NACA 0012 of xfoil-python/xfoil/test.py, plus the recorded baseline values (API_SNAPSHOT.md)."""
from ctypes import byref, c_bool, c_double, c_int, c_void_p

import numpy as np
import pytest

from conftest import same_reals, sync_globals
from test_xbl import diff_allowing_nan_states



def naca0012():
    """x, y arrays of the Airfoil(...) literal in the original test file (embedded copy)."""
    from golden import naca0012_coords as c
    return np.array(c.X, dtype=np.float64), np.array(c.Y, dtype=np.float64)


def ptr(a):
    return a.ctypes.data_as(c_void_p)


def setup(clib, flib, re=1e6, mach=0.0, max_iter=100, repanel=None):
    sync_globals(clib, flib)
    x, y = naca0012()
    for lib in (clib, flib):
        lib.init()
        lib.set_airfoil(ptr(x.copy()), ptr(y.copy()), byref(c_int(len(x))))
        if repanel:
            lib.repanel(byref(c_int(repanel)), *[byref(c_double(v)) for v in (1.0, 0.15, 0.2, 1, 1, 1, 1)])
        lib.set_reynolds(byref(c_double(re)))
        lib.set_mach(byref(c_double(mach)))
        lib.set_max_iter(byref(c_int(max_iter)))
    assert diff_allowing_nan_states(clib, flib) == []


def point(lib, name, value):
    out = [c_double() for _ in range(6)]
    conv = c_bool()
    getattr(lib, name)(byref(c_double(value)), *[byref(o) for o in out], byref(conv))
    return np.array([o.value for o in out]), conv.value


def seq(lib, name, start, end, step):
    n = abs(int((end - start) / step))                      # as xfoil.py
    arrs = [np.zeros(n) for _ in range(5)]
    conv = np.zeros(n, dtype=np.bool_)
    getattr(lib, name)(byref(c_double(start)), byref(c_double(end)), byref(c_int(n)),
                       *[ptr(a) for a in arrs], ptr(conv))
    return np.concatenate(arrs), conv


def cp_dist(lib):
    lib.get_n_cp.restype = c_int
    n = lib.get_n_cp()
    arrs = [np.zeros(n) for _ in range(3)]
    lib.get_cp(*[ptr(a) for a in arrs], byref(c_int(n)))
    return np.concatenate(arrs)


def compare_point(clib, flib, name, value):
    (c, cc), (f, fc) = point(clib, name, value), point(flib, name, value)
    assert same_reals(c, f) and cc == fc, (name, value, c, f)
    assert diff_allowing_nan_states(clib, flib) == []
    cd_, fd_ = cp_dist(clib), cp_dist(flib)
    n = len(cd_) // 3
    assert same_reals(cd_[:2 * n], fd_[:2 * n])                 # x, y
    if cc:
        assert same_reals(cd_[2 * n:], fd_[2 * n:])
    else:                                                       # clean-up C6c: no cp of an unconverged point
        assert np.isnan(cd_[2 * n:]).all()
    return c, cc


def test_a_baseline(clib, flib):
    setup(clib, flib, max_iter=100)
    (cl, cd, cm, cp, _, _), conv = compare_point(clib, flib, 'alfa', 10.0)
    assert conv
    assert np.round([cl, cd, cm, cp], 4).tolist() == [1.0809, 0.0150, 0.0053, -5.5077]


def test_cl_baseline(clib, flib):
    setup(clib, flib, max_iter=40)
    (a, cd, cm, cp, _, _), conv = compare_point(clib, flib, 'cl', 1.0)
    assert conv
    assert np.round([a, cd, cm, cp], 4).tolist() == [9.0617, 0.0135, 0.0013, -4.7361]


@pytest.mark.parametrize('name,start,end,step,max_iter', [('aseq', -20, 20, 0.5, 150),
                                                          ('cseq', -0.5, 0.5, 0.05, 150)])
def test_sequences_match_fortran(clib, flib, name, start, end, step, max_iter):
    setup(clib, flib, max_iter=max_iter)
    (c, cc), (f, fc) = seq(clib, name, start, end, step), seq(flib, name, start, end, step)
    assert same_reals(c, f) and np.array_equal(cc, fc)
    assert diff_allowing_nan_states(clib, flib) == []
    assert cc.any()


@pytest.mark.parametrize('re,mach,panels', [(0.0, 0.0, None), (0.0, 0.4, 120), (3e5, 0.3, 200), (2e6, 0.0, 90)])
def test_points_and_short_sequences_match_fortran(clib, flib, re, mach, panels):
    setup(clib, flib, re=re, mach=mach, max_iter=60, repanel=panels)
    for name, v in (('alfa', 0.0), ('alfa', 4.5), ('cl', 0.6), ('alfa', -3.0), ('cl', 1.4)):
        compare_point(clib, flib, name, v)
    for name, args in (('aseq', (-4, 8, 2.0)), ('cseq', (0.2, 1.0, 0.2))):
        (c, cc), (f, fc) = seq(clib, name, *args), seq(flib, name, *args)
        assert same_reals(c, f) and np.array_equal(cc, fc), name
        assert diff_allowing_nan_states(clib, flib) == []


@pytest.mark.parametrize('name,args', [('aseq', (-6, 18, 1.0)), ('cseq', (0.0, 1.6, 0.1))])
def test_unconverged_points_match_fortran(clib, flib, name, args):
    """Few iterations: some points do not converge (conv False -> NaN in Python)."""
    setup(clib, flib, re=5e5, max_iter=1)
    (c, cc), (f, fc) = seq(clib, name, *args), seq(flib, name, *args)
    assert same_reals(c, f) and np.array_equal(cc, fc)
    assert diff_allowing_nan_states(clib, flib) == []
    assert not cc.all(), cc
    compare_point(clib, flib, 'alfa', 14.0)
