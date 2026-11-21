# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""m_xbl.f90 -> src/xbl.c. Part 1: dslim, blpini, xifset (side argument: C 0-based, Fortran 1-based)."""
import os
import subprocess
import sys
from ctypes import byref, c_double, c_int, c_void_p

import numpy as np
import pytest

from conftest import global_array, globals_diff, sync_globals, va_vb_diff

IVX, IZX = 452, 722


def test_dslim_matches_fortran(clib, flib):
    rng = np.random.default_rng(0)
    for _ in range(3000):
        th = 10.0 ** rng.uniform(-5, -2)
        args = (float(th * rng.uniform(1.0, 6.0)), float(th), float(rng.uniform(0.1, 2)),
                float(rng.uniform(0, 0.8)), float(rng.uniform(1.02, 2.5)))
        cd, fd = c_double(args[0]), c_double(args[0])
        clib.dslim(byref(cd), *[c_double(a) for a in args[1:]])
        flib.__m_xbl_MOD_dslim(byref(fd), *[byref(c_double(a)) for a in args[1:]])
        assert cd.value.hex() == fd.value.hex(), args


def test_blpini_matches_fortran(clib, flib):
    sync_globals(clib, flib)
    for lib, fortran in ((clib, False), (flib, True)):
        global_array(lib, 'i_blpar:ctcon', np.float64, None, fortran)[:] = -1.0
    clib.blpini()
    flib.__m_xbl_MOD_blpini()
    assert globals_diff(clib, flib) == []
    assert global_array(clib, 'i_blpar:gacon', np.float64, None, False)[0] == 6.70


def setup_xifset(clib, flib, seed):
    import test_xpanel as P
    st, rng, n, nw = P.setup_state(clib, flib, seed)
    x, _ = st.views('x'); y, _ = st.views('y'); s, _ = st.views('s')
    ile = int(np.argmin(x[:n]))
    st.set('sle', float(s[ile])); st.set('xle', float(x[ile])); st.set('yle', float(y[ile]))
    st.set('xte', 0.5 * (x[0] + x[n - 1])); st.set('yte', 0.5 * (y[0] + y[n - 1]))
    chord = float(np.ptp(x[:n]))
    st.set('sst', float(s[ile]) + float(rng.choice([0.0, 0.002, -0.002, 0.3, -0.3])) * chord)
    st.set_index('iblte', [ile + 2, n - ile])
    xssi = np.zeros((IVX, 2), order='F')
    xssi[:, :] = rng.uniform(0.5, 2.0, (IVX, 2)) * chord
    st.set('xssi', xssi)
    return st, rng


def test_xifset_matches_fortran(clib, flib):
    sync_globals(clib, flib)
    hits = {'te': 0, 'past': 0, 'normal': 0}
    for seed in range(24):
        st, rng = setup_xifset(clib, flib, seed)
        for is_ in (0, 1):
            for xs in (0.0, 0.05, 0.3, 0.9, 1.0, 1.5):
                st.set('xstrip', [xs, xs])
                clib.xifset(c_int(is_))
                flib.__m_xbl_MOD_xifset(byref(c_int(is_ + 1)))
                assert globals_diff(clib, flib) == [], (seed, is_, xs)
                xif = global_array(clib, 'i_xbl:xiforc', np.float64, None, False)[0]
                xte_bl = st.views('xssi')[0].reshape((IVX, 2), order='F')[st.views('iblte', np.int32)[0][is_], is_]
                hits['te' if xs >= 1.0 else 'past' if xif == xte_bl else 'normal'] += 1
    assert min(hits.values()) > 10, hits


# ------------------------------------------------------------------ part 2: mrchue

def prepare_march(clib, flib, seed, closed_te=False, sharp_le=False):
    """Panel + wake + BL-mapping chain (verified in test_xpanel), inviscid Ue, BL parameters."""
    import test_xblsys as B
    import test_xpanel as P
    st, n = P.run_wake_chain(clib, flib, seed, closed_te, sharp_le)
    for step in ('qwcalc', 'qiset', 'uicalc'):
        getattr(clib, step)()
        getattr(flib, '__m_xpanel_MOD_' + step)()
    rng = np.random.default_rng(30000 + seed)
    B.set_bl_parameters(clib, flib, rng)
    st.set('uedg', st.views('uinv')[0].copy())
    st.set('acrit', [float(rng.uniform(5, 11))] * 2)
    st.set('xstrip', [float(rng.choice([1.0, rng.uniform(0.05, 0.9)])) for _ in range(2)])
    x, _ = st.views('x'); y, _ = st.views('y'); s, _ = st.views('s')
    ile = int(np.argmin(x[:n]))
    st.set('sle', float(s[ile])); st.set('xle', float(x[ile])); st.set('yle', float(y[ile]))
    idampv = int(rng.integers(0, 2))
    for lib, fortran in ((clib, False), (flib, True)):
        global_array(lib, 'i_xbl:idampv', np.int32, None, fortran)[:] = idampv
    assert globals_diff(clib, flib) == []
    return st, n


import test_xpanel as _P  # noqa: E402

MARCH_STATES = _P.STATES


@pytest.mark.parametrize('seed,closed_te,sharp_le', MARCH_STATES)
def test_mrchue_matches_fortran(clib, flib, seed, closed_te, sharp_le):
    sync_globals(clib, flib)
    prepare_march(clib, flib, seed, closed_te, sharp_le)
    clib.mrchue.restype = c_int
    flib.__m_xbl_MOD_mrchue.restype = c_int
    rc = clib.mrchue()
    rf = flib.__m_xbl_MOD_mrchue()
    assert rc == rf
    assert globals_diff(clib, flib) == []


# ------------------------------------------------------------------ part 3: mrchdu

@pytest.mark.parametrize('seed,closed_te,sharp_le', MARCH_STATES)
def test_mrchdu_matches_fortran(clib, flib, seed, closed_te, sharp_le):
    """mrchue first (initial BL, as XFOIL's setbl does), then mrchdu from that state, twice."""
    sync_globals(clib, flib)
    prepare_march(clib, flib, seed, closed_te, sharp_le)
    for lib, name in ((clib, 'mrchue'), (flib, '__m_xbl_MOD_mrchue'), (clib, 'mrchdu'), (flib, '__m_xbl_MOD_mrchdu')):
        getattr(lib, name).restype = c_int
    assert clib.mrchue() == flib.__m_xbl_MOD_mrchue()
    assert globals_diff(clib, flib) == []
    for rep in range(2):
        rng = np.random.default_rng(seed + 100 * rep)
        uc = global_array(clib, 'uedg', np.float64, (IVX, 2), False)
        uf = global_array(flib, 'uedg', np.float64, (IVX, 2), True)
        bump = 1.0 + 0.02 * rng.normal(size=uc.shape)            # perturb Ue as a Newton update would
        uc[...] = uc * bump
        uf[...] = uf * bump
        assert clib.mrchdu() == flib.__m_xbl_MOD_mrchdu(), rep
        assert globals_diff(clib, flib) == [], rep


# ------------------------------------------------------------------ part 4: setbl (+ blsolv: one Newton step)

def prepare_setbl(clib, flib, seed, closed_te=False, sharp_le=False):
    """prepare_march + source influence (qdcalc), BL system lines (iblsys), Mach/Re(CL) setup."""
    st, n = prepare_march(clib, flib, seed, closed_te, sharp_le)
    clib.qdcalc(); flib.__m_xpanel_MOD_qdcalc()
    clib.iblsys(); flib.__s_xbl_MOD_iblsys()
    rng = np.random.default_rng(40000 + seed)
    reybl = global_array(clib, 'i_xbl:reybl', np.float64, None, False)[0]
    vals = {'gamma': 1.4, 'gamm1': 0.4, 'minf1': float(rng.choice([0.0, rng.uniform(0.05, 0.6)])),
            'reinf1': float(reybl), 'cl': float(rng.uniform(-0.5, 1.5)), 'clspec': float(rng.uniform(0.1, 1.2))}
    for name, v in vals.items():
        st.set(name, v)
    for name, v in (('matyp', 1), ('retyp', 1), ('lalfa', int(seed % 2 == 0)), ('lblini', 0),
                    ('idamp', int(rng.integers(0, 2)))):
        st.set(name, v, np.int32)
    assert globals_diff(clib, flib) == []
    return st, n


def diff_allowing_nan_states(clib, flib):
    """globals_diff, but in diverged states (NaN in thet) only require identical NaN positions
    of the BL arrays (user priority: no chasing of compiler NaN behaviour)."""
    d = globals_diff(clib, flib)
    th_c = global_array(clib, 'thet', np.float64, (IVX, 2), False)
    if not d or not np.isnan(th_c).any():
        return d
    for name in ('thet', 'dstr', 'uedg', 'mass'):
        c = global_array(clib, name, np.float64, (IVX, 2), False)
        f = global_array(flib, name, np.float64, (IVX, 2), True)
        assert np.array_equal(np.isnan(c), np.isnan(f)), name
    return []


@pytest.mark.parametrize('seed,closed_te,sharp_le', MARCH_STATES)
def test_setbl_and_newton_step_match_fortran(clib, flib, seed, closed_te, sharp_le):
    sync_globals(clib, flib)
    prepare_setbl(clib, flib, seed, closed_te, sharp_le)
    clib.setbl.restype = c_int
    flib.__m_xbl_MOD_setbl.restype = c_int
    for rep in range(2):                        # 1st: mrchue + mrchdu; 2nd: lblini set, mrchdu only
        assert clib.setbl() == flib.__m_xbl_MOD_setbl(), rep
        assert diff_allowing_nan_states(clib, flib) == [], ('setbl', rep)
        assert va_vb_diff(clib, flib) == [], ('setbl', rep)
    clib.blsolv(); flib.__m_xsolve_MOD_blsolv()
    assert diff_allowing_nan_states(clib, flib) == [], 'blsolv'
    assert va_vb_diff(clib, flib) == [], 'blsolv'


# ------------------------------------------------------------------ part 5: update (full Newton iterations)

@pytest.mark.parametrize('seed,closed_te,sharp_le', MARCH_STATES)
def test_newton_iterations_with_update_match_fortran(clib, flib, seed, closed_te, sharp_le):
    st, _ = prepare_setbl(clib, flib, seed, closed_te, sharp_le)
    rng = np.random.default_rng(50000 + seed)
    st.set('dtor', np.pi / 180.0)
    st.set('minf_cl', float(rng.uniform(-0.2, 0.2)))
    st.set('matyp', int(1 + seed % 3 // 2), np.int32)       # matyp 1 and 2 (dclmin branch)
    assert globals_diff(clib, flib) == []
    clib.setbl.restype = c_int
    flib.__m_xbl_MOD_setbl.restype = c_int
    for it in range(3):
        ok = clib.setbl()
        assert ok == flib.__m_xbl_MOD_setbl(), it
        if ok:                                  # aborted (NaN) setbl: VA/VB not rebuilt, hold leftovers
            assert va_vb_diff(clib, flib) == [], ('setbl', it)
        clib.blsolv(); flib.__m_xsolve_MOD_blsolv()
        clib.update(); flib.__m_xbl_MOD_update()
        f_rlx = float(global_array(flib, 'rlx', np.float64, None, True)[0])
        if not 0.0 < f_rlx <= 1.0:
            # clean-up C6c: Fortran's DCLMIN limit gave a negative/huge RLX (Re(CL) and CL < 0); the C
            # limit only reduces RLX, so the states part here
            c_rlx = float(global_array(clib, 'rlx', np.float64, None, False)[0])
            assert 0.0 < c_rlx <= 1.0, c_rlx
            break
        assert diff_allowing_nan_states(clib, flib) == [], ('update', it)
        st.set('lblini', 1, np.int32)
    vmx = global_array(clib, 'vmxbl', np.uint8, None, False)[0]
    assert chr(vmx) in 'nCTDU '                 # blank: no finite change (diverged state)
