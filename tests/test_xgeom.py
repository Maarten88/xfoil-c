# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""m_xgeom.f90 (API-reachable part) -> src/xgeom.c: results must be bit-identical to Fortran.

Test shapes are NACA airfoils from the (already verified) C naca4/naca5, with random
point counts, rotation, scaling, translation, noise and optional sharp (doubled) LE point.
"""
import subprocess
import sys
from ctypes import byref, c_double, c_int, c_void_p, create_string_buffer

import numpy as np
import pytest

from conftest import FORTRAN_LIB

N_SHAPES = 150


def ptr(arr):
    return arr.ctypes.data_as(c_void_p)


def F(flib, name):
    return getattr(flib, '__m_xgeom_MOD_' + name)


def naca_shape(clib, rng, sharp_le=False):
    nside = int(rng.integers(20, 214))
    xx, yt, yc = np.zeros(nside), np.zeros(nside), np.zeros(nside)
    xb, yb, nb = np.zeros(2 * nside), np.zeros(2 * nside), c_int()
    name = create_string_buffer(48)
    if rng.random() < 0.7:
        ides = int(rng.integers(0, 10)) * 1000 + int(rng.integers(0, 10)) * 100 + int(rng.integers(6, 25))
        clib.naca4(c_int(ides), ptr(xx), ptr(yt), ptr(yc), c_int(nside), ptr(xb), ptr(yb), byref(nb),
                   name, c_int(48))
    else:
        ides = c_int(int(rng.choice([210, 220, 230, 240, 250])) * 100 + int(rng.integers(6, 25)))
        clib.naca5(byref(ides), ptr(xx), ptr(yt), ptr(yc), c_int(nside), ptr(xb), ptr(yb), byref(nb),
                   name, c_int(48))
    x, y = xb[:nb.value].copy(), yb[:nb.value].copy()
    if sharp_le:                                  # doubled LE point -> segment joint in s
        x, y = np.insert(x, nside - 1, x[nside - 1]), np.insert(y, nside - 1, y[nside - 1])
    if rng.random() < 0.7:                        # rotate, scale, translate, add noise
        a = float(rng.uniform(-0.5, 0.5))
        c, s = np.cos(a), np.sin(a)
        x, y = c * x - s * y, s * x + c * y
        k = float(rng.uniform(0.1, 10.0))
        x, y = k * x + float(rng.uniform(-5, 5)), k * y + float(rng.uniform(-5, 5))
        if not sharp_le:
            x, y = x + rng.normal(0, 1e-5 * k, x.size), y + rng.normal(0, 1e-5 * k, y.size)
    return np.ascontiguousarray(x), np.ascontiguousarray(y)


def wavy_shape(clib, seed):
    """Noisy, wavy NACA shape: makes lefind/sopps hit their failure fallbacks
    (e.g. seeds 9, 28, 44, 99 for lefind; 9, 16, 20, 26 for sopps)."""
    rng = np.random.default_rng(seed)
    x, y = naca_shape(clib, rng)
    amp = float(rng.choice([1e-3, 3e-3, 1e-2, 3e-2]))
    k = float(rng.uniform(5, 60))
    t = np.linspace(0, 1, len(x))
    y = y + amp * np.sin(k * 2 * np.pi * t) + rng.normal(0, amp / 3, len(x))
    return np.ascontiguousarray(x), np.ascontiguousarray(y)


WAVY_SEEDS = sorted(set(range(100)) | {99, 162, 261, 417, 432, 451})


def splined(clib, x, y):
    n = len(x)
    s, xp, yp = np.zeros(n), np.zeros(n), np.zeros(n)
    clib.scalc(ptr(x), ptr(y), ptr(s), c_int(n))
    clib.segspl(ptr(x), ptr(xp), ptr(s), c_int(n))
    clib.segspl(ptr(y), ptr(yp), ptr(s), c_int(n))
    return s, xp, yp


def shapes(clib, n=N_SHAPES, seed0=0):
    for seed in range(seed0, seed0 + n):
        rng = np.random.default_rng(seed)
        x, y = naca_shape(clib, rng, sharp_le=(seed % 5 == 0))
        s, xp, yp = splined(clib, x, y)
        yield rng, x, y, s, xp, yp


def bits(*vals):
    return tuple(float(v).hex() for v in vals)


def run_lefind(lib, fortran, x, xp, y, yp, s):
    sle = c_double(0.0)
    if fortran:
        F(lib, 'lefind')(byref(sle), ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s), byref(c_int(len(x))))
    else:
        lib.lefind(byref(sle), ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s), c_int(len(x)))
    return sle.value


def test_lefind_matches_fortran(clib, flib):
    sharp = 0
    for _, x, y, s, xp, yp in shapes(clib):
        c = run_lefind(clib, False, x, xp, y, yp, s)
        assert bits(c) == bits(run_lefind(flib, True, x, xp, y, yp, s))
        sharp += bool(np.any(s[c == s][:1]) and np.sum(s == c) > 1)
    assert sharp > 0                                  # sharp-LE early return exercised


def test_sopps_matches_fortran(clib, flib):
    for _, x, y, s, xp, yp in shapes(clib, n=30):
        n = len(x)
        sle = run_lefind(clib, False, x, xp, y, yp, s)
        for si in np.concatenate([s, [sle]]):
            cs, fs = c_double(), c_double()
            clib.sopps(byref(cs), c_double(si), ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s), c_int(n),
                       c_double(sle))
            F(flib, 'sopps')(byref(fs), byref(c_double(si)), ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s),
                             byref(c_int(n)), byref(c_double(sle)))
            assert bits(cs.value) == bits(fs.value)


def test_norm_matches_fortran(clib, flib):
    for _, x, y, _, _, _ in shapes(clib):
        n = len(x)
        out = []
        for fortran in (False, True):
            ax, ay = x.copy(), y.copy()
            s, xp, yp = np.zeros(n), np.zeros(n), np.zeros(n)
            if fortran:
                F(flib, 'norm')(ptr(ax), ptr(xp), ptr(ay), ptr(yp), ptr(s), byref(c_int(n)))
            else:
                clib.norm(ptr(ax), ptr(xp), ptr(ay), ptr(yp), ptr(s), c_int(n))
            out.append(b''.join(a.tobytes() for a in (ax, ay, s, xp, yp)))
        assert out[0] == out[1]


def run_geopar(lib, fortran, x, xp, y, yp, s):
    n = len(x)
    t = np.zeros(n)
    outs = [c_double(-1.0) for _ in range(15)]
    if fortran:
        F(lib, 'geopar')(ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s), byref(c_int(n)), ptr(t),
                         *[byref(o) for o in outs])
    else:
        lib.geopar(ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s), c_int(n), ptr(t), *[byref(o) for o in outs])
    return bits(*[o.value for o in outs]), t.tobytes()


def test_geopar_matches_fortran(clib, flib):
    for _, x, y, s, xp, yp in shapes(clib):
        assert run_geopar(clib, False, x, xp, y, yp, s) == run_geopar(flib, True, x, xp, y, yp, s)


def run_aecalc(lib, fortran, x, y, t, itype):
    n = len(x)
    outs = [c_double(-1.0) for _ in range(7)]
    if fortran:
        F(lib, 'aecalc')(byref(c_int(n)), ptr(x), ptr(y), ptr(t), byref(c_int(itype)), *[byref(o) for o in outs])
    else:
        lib.aecalc(c_int(n), ptr(x), ptr(y), ptr(t), c_int(itype), *[byref(o) for o in outs])
    return bits(*[o.value for o in outs])


@pytest.mark.parametrize('itype', [1, 2])
def test_aecalc_matches_fortran(clib, flib, itype):
    for rng, x, y, _, _, _ in shapes(clib):
        t = rng.uniform(0.001, 0.01, len(x))
        assert run_aecalc(clib, False, x, y, t, itype) == run_aecalc(flib, True, x, y, t, itype)


@pytest.mark.parametrize('case', ['collinear', 'zero_t', 'circle', 'square'])
def test_aecalc_special_sections(clib, flib, case):
    th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    x, y, t = np.cos(th), np.sin(th), np.full(64, 0.01)
    if case == 'collinear':
        x, y = np.linspace(0, 1, 64), np.zeros(64)
    elif case == 'zero_t':
        t = np.zeros(64)
    elif case == 'square':
        x, y = np.array([0.0, 1.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0, 1.0])
        t = np.full(4, 0.01)
    for itype in (1, 2):
        assert run_aecalc(clib, False, x, y, t, itype) == run_aecalc(flib, True, x, y, t, itype)


def test_tccalc_matches_fortran(clib, flib):
    for _, x, y, s, xp, yp in shapes(clib):
        n = len(x)
        res = []
        for fortran in (False, True):
            o = [c_double() for _ in range(4)]
            if fortran:
                F(flib, 'tccalc')(ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s), byref(c_int(n)), *[byref(v) for v in o])
            else:
                clib.tccalc(ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s), c_int(n), *[byref(v) for v in o])
            res.append(bits(*[v.value for v in o]))
        assert res[0] == res[1]


def run_cang(lib, fortran, x, y, iprint=0):
    n = len(x)
    imax, amax = c_int(-5), c_double()
    if fortran:
        F(lib, 'cang')(ptr(x), ptr(y), byref(c_int(n)), byref(c_int(iprint)), byref(imax), byref(amax))
        return imax.value - 1, bits(amax.value)          # Fortran index is 1-based
    lib.cang(ptr(x), ptr(y), c_int(n), c_int(iprint), byref(imax), byref(amax))
    return imax.value, bits(amax.value)


def test_cang_matches_fortran(clib, flib):
    doubled = 0
    for rng, x, y, _, _, _ in shapes(clib):
        if rng.random() < 0.5 and len(x) > 10:          # interior doubled point
            k = int(rng.integers(3, len(x) - 3))
            x, y = np.insert(x, k, x[k]), np.insert(y, k, y[k])
            doubled += 1
        assert run_cang(clib, False, x, y) == run_cang(flib, True, x, y)
    assert doubled > 20


def test_failure_fallbacks_match_fortran(clib, flib):
    """lefind (50 Newton steps) and sopps (12 steps / resd == 0) fallbacks on wavy shapes."""
    for seed in WAVY_SEEDS:
        x, y = wavy_shape(clib, seed)
        s, xp, yp = splined(clib, x, y)
        n = len(x)
        sle = run_lefind(clib, False, x, xp, y, yp, s)
        assert bits(sle) == bits(run_lefind(flib, True, x, xp, y, yp, s)), seed
        for si in s:
            cs, fs = c_double(), c_double()
            clib.sopps(byref(cs), c_double(si), ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s), c_int(n),
                       c_double(sle))
            F(flib, 'sopps')(byref(fs), byref(c_double(si)), ptr(x), ptr(xp), ptr(y), ptr(yp), ptr(s),
                             byref(c_int(n)), byref(c_double(sle)))
            assert bits(cs.value) == bits(fs.value), (seed, si)
        assert run_geopar(clib, False, x, xp, y, yp, s) == run_geopar(flib, True, x, xp, y, yp, s), seed


# ------------------------------------------------------------------ printed output


