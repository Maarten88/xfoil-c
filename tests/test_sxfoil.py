# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""s_xfoil.f90 -> src/sxfoil.c: outputs and the whole global state identical to Fortran."""
import math
import os
import subprocess
import sys
from ctypes import byref, c_double, c_int, c_void_p

import numpy as np
import pytest

from conftest import global_array, globals_diff, same_reals, sync_globals


def ptr(arr):
    return arr.ctypes.data_as(c_void_p)


def set_both(clib, flib, name, value, dtype=np.float64):
    for lib, fortran in ((clib, False), (flib, True)):
        global_array(lib, name, dtype, None, fortran)[:] = value


MRCL_CASES = [(m, r, cls, minf1, reinf1)
              for m in (0, 1, 2, 3, 4) for r in (0, 1, 2, 3, 4)
              for cls in (1.0, 0.3, 1e-5, -0.2)
              for minf1, reinf1 in ((0.0, 1e6), (0.5, 3e5), (0.3, 0.0))]


def run_mrcl(clib, flib, case):
    m, r, cls, minf1, reinf1 = case
    set_both(clib, flib, 'matyp', m, np.int32)
    set_both(clib, flib, 'retyp', r, np.int32)
    set_both(clib, flib, 'minf1', minf1)
    set_both(clib, flib, 'reinf1', reinf1)
    cm, cr, fm, fr = c_double(), c_double(), c_double(), c_double()
    clib.mrcl(c_double(cls), byref(cm), byref(cr))
    flib.__s_xfoil_MOD_mrcl(byref(c_double(cls)), byref(fm), byref(fr))
    return (cm.value, cr.value), (fm.value, fr.value)


def test_mrcl_matches_fortran(clib, flib):
    sync_globals(clib, flib)
    for case in MRCL_CASES:
        c, f = run_mrcl(clib, flib, case)
        assert same_reals(np.array(c), np.array(f)), case
        assert globals_diff(clib, flib) == [], case


def test_comset_matches_fortran(clib, flib):
    sync_globals(clib, flib)
    set_both(clib, flib, 'gamma', 1.4)
    set_both(clib, flib, 'gamm1', 0.4)
    for minf in [0.0, 0.05, 0.3, 0.7, 0.95] + list(np.random.default_rng(0).uniform(0, 0.99, 30)):
        set_both(clib, flib, 'minf', float(minf))
        set_both(clib, flib, 'qinf', float(1.0 + minf))
        clib.comset()
        flib.__s_xfoil_MOD_comset()
        assert globals_diff(clib, flib) == [], minf


def test_cpcalc_matches_fortran(clib, flib):
    for seed in range(100):
        rng = np.random.default_rng(seed)
        n = int(rng.integers(1, 300))
        q = rng.uniform(-3.0, 3.0, n)
        qinf, minf = float(rng.uniform(0.5, 2.0)), float(rng.choice([0.0, rng.uniform(0, 0.9)]))
        cc, cf = np.zeros(n), np.zeros(n)
        clib.cpcalc(c_int(n), ptr(q), c_double(qinf), c_double(minf), ptr(cc))
        flib.__s_xfoil_MOD_cpcalc(byref(c_int(n)), ptr(q), byref(c_double(qinf)), byref(c_double(minf)), ptr(cf))
        assert same_reals(cc, cf), seed


def test_clcalc_and_tecalc_match_fortran(clib, flib):
    import test_xpanel as P
    for seed, closed in [(s, c) for s in range(8) for c in (False, True)]:
        sync_globals(clib, flib)
        st, rng, n, nw = P.setup_state(clib, flib, seed, closed_te=closed)
        x, _ = st.views('x'); y, _ = st.views('y')
        gam = rng.normal(size=n); gam_a = rng.normal(size=n)
        st.set('gam', gam); st.set('gam_a', gam_a)
        st.set('chord', float(np.ptp(x[:n])))
        clib.tecalc(); flib.__s_xfoil_MOD_tecalc()
        assert globals_diff(clib, flib) == [], ('tecalc', seed, closed)
        xs, ys = np.ascontiguousarray(x[:n]), np.ascontiguousarray(y[:n])
        for alfa, minf in ((0.0, 0.0), (0.1, 0.3), (-0.2, 0.6)):
            args = (xs, ys, gam, gam_a)
            sc = (alfa, minf, 1.2, 0.25, 0.01)
            co = [c_double() for _ in range(5)]
            fo = [c_double() for _ in range(5)]
            clib.clcalc(c_int(n), *[ptr(a) for a in args], *[c_double(v) for v in sc], *[byref(o) for o in co])
            flib.__s_xfoil_MOD_clcalc(byref(c_int(n)), *[ptr(a) for a in args], *[byref(c_double(v)) for v in sc],
                                      *[byref(o) for o in fo])
            assert same_reals(np.array([o.value for o in co]), np.array([o.value for o in fo])), (seed, alfa, minf)


def test_cdcalc_matches_fortran(clib, flib):
    import test_xbl as XB
    for seed in range(6):
        sync_globals(clib, flib)
        XB.prepare_march(clib, flib, seed)
        clib.mrchue(); flib.__m_xbl_MOD_mrchue()
        for lvisc, lblini in ((1, 1), (0, 1), (1, 0)):
            set_both(clib, flib, 'lvisc', lvisc, np.int32)
            set_both(clib, flib, 'lblini', lblini, np.int32)
            set_both(clib, flib, 'tklam', float(np.random.default_rng(seed).uniform(0, 0.1)))
            clib.cdcalc(); flib.__s_xfoil_MOD_cdcalc()
            assert globals_diff(clib, flib) == [], (seed, lvisc, lblini)


