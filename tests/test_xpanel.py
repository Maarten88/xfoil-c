# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""m_xpanel.f90 (apcalc, ncalc, psilin, pswlin) -> src/xpanel.c.

Identical i_xfoil state is written into both libraries (airfoil from naca_shape, a wake,
random vorticity/source strengths); after each call the outputs and every global the
routine writes must be byte-identical. Node indices: C 0-based (-1 = off-body), Fortran
1-based (0 = off-body).
"""
import math
from ctypes import byref, c_char, c_double, c_int, c_size_t, c_void_p

import numpy as np
import pytest

import test_xgeom as G

IQX, IZX = 640, 722

WRITTEN = ['dzdg', 'dzdn', 'dqdg', 'dzdm', 'dqdm', 'z_qinf', 'z_alfa', 'z_qdof0', 'z_qdof1',
           'z_qdof2', 'z_qdof3', 'qtan1', 'qtan2', 'sigte', 'gamte', 'cosa', 'sina', 'apanel']


def ptr(arr):
    return arr.ctypes.data_as(c_void_p)


class State:
    """Byte views of the same i_xfoil globals in the C and the Fortran library."""

    def __init__(self, clib, flib):
        self.clib, self.flib = clib, flib
        clib.xfoil_global.restype = c_void_p

    def views(self, name, dtype=np.float64):
        size = c_size_t()
        addr = self.clib.xfoil_global(name.encode(), byref(size))
        c = np.frombuffer((c_char * size.value).from_address(addr), dtype=dtype)
        f = np.frombuffer((c_char * size.value).in_dll(self.flib, '__i_xfoil_MOD_' + name), dtype=dtype)
        return c, f

    def set(self, name, value, dtype=np.float64):
        for v in self.views(name, dtype):
            v[:] = 0
            val = np.asarray(value, dtype=dtype).ravel(order='F')
            v[:val.size] = val

    def set2(self, name, cols, ld, dtype=np.float64):
        """2-D global (ld rows in Fortran order): cols = list of 1-D columns."""
        arr = np.zeros((ld, len(cols)), dtype=dtype, order='F')
        for k, col in enumerate(cols):
            arr[:len(col), k] = col
        self.set(name, arr, dtype)

    def set_index(self, name, fortran_values):
        """Index global: Fortran gets the 1-based values, C the mapped 0-based ones."""
        from conftest import to_c_index
        f_vals = np.asarray(fortran_values, dtype=np.int32)
        c, f = self.views(name, np.int32)
        f[:f_vals.size] = f_vals
        c[:f_vals.size] = to_c_index(f_vals)

    def same(self, name, dtype=np.float64):
        from conftest import same_reals
        c, f = self.views(name, dtype)
        return same_reals(c, f) if dtype == np.float64 else c.tobytes() == f.tobytes()


def setup_state(clib, flib, seed, closed_te=False, sharp_le=False):
    rng = np.random.default_rng(seed)
    st = State(clib, flib)
    x, y = G.naca_shape(clib, rng, sharp_le=sharp_le)
    n = len(x)
    if closed_te:
        x[-1], y[-1] = x[0], y[0]
    s, xp, yp = G.splined(clib, x, y)
    xn, yn = np.zeros(n), np.zeros(n)
    clib.ncalc(ptr(x), ptr(y), ptr(s), c_int(n), ptr(xn), ptr(yn))
    nw = int(rng.integers(5, 40))
    xte, yte = 0.5 * (x[0] + x[-1]), 0.5 * (y[0] + y[-1])
    chord = float(np.ptp(x))
    wx = xte + np.cumsum(rng.uniform(0.01, 0.08, nw)) * chord
    wy = yte + (wx - xte) * float(rng.uniform(-0.1, 0.1))
    ang = rng.uniform(0, 2 * np.pi, nw)
    pi = 4.0 * math.atan(1.0)
    st.set('n', n, np.int32)
    st.set('nw', nw, np.int32)
    st.set('x', np.concatenate([x, wx]))
    st.set('y', np.concatenate([y, wy]))
    st.set('s', np.concatenate([s, s[-1] + np.cumsum(np.hypot(np.diff(np.r_[xte, wx]), np.diff(np.r_[yte, wy])))]))
    st.set('xp', xp)
    st.set('yp', yp)
    st.set('xte', xte)
    st.set('yte', yte)
    st.set('nx', np.concatenate([xn, np.cos(ang)]))
    st.set('ny', np.concatenate([yn, np.sin(ang)]))
    st.set('pi', pi)
    st.set('qopi', 0.25 / pi)
    st.set('hopi', 0.50 / pi)
    st.set('qinf', float(rng.uniform(0.5, 2.0)))
    st.set('alfa', float(rng.uniform(-0.3, 0.3)))
    st.set('sharp', int(closed_te), np.int32)
    st.set('ante', float(rng.uniform(1e-3, 1e-2)))
    st.set('aste', float(rng.uniform(-1e-3, 1e-3)))
    st.set('dste', float(rng.uniform(2e-3, 2e-2)))
    st.set('yimage', float(rng.uniform(-2.0, -0.5)))
    st.set('gam', rng.normal(size=n))
    st.set2('gamu', [rng.normal(size=n), rng.normal(size=n)], IQX)
    st.set('sig', rng.normal(size=n + nw) * 0.1)
    for q in ('qf0', 'qf1', 'qf2', 'qf3'):
        st.set(q, rng.normal(size=n))
    for name in WRITTEN:                       # identical stale contents
        st.set(name, rng.normal(size=st.views(name)[0].size))
    return st, rng, n, nw


def call_psilin(lib, fortran, i, xi, yi, nxi, nyi, geolin, siglin):
    psi, psi_ni = c_double(), c_double()
    if fortran:
        lib.__m_xpanel_MOD_psilin(byref(c_int(i + 1)), byref(c_double(xi)), byref(c_double(yi)),
                                  byref(c_double(nxi)), byref(c_double(nyi)), byref(psi), byref(psi_ni),
                                  byref(c_int(geolin)), byref(c_int(siglin)))
    else:
        lib.psilin(c_int(i), c_double(xi), c_double(yi), c_double(nxi), c_double(nyi), byref(psi),
                   byref(psi_ni), c_int(geolin), c_int(siglin))
    return psi.value.hex(), psi_ni.value.hex()


def call_pswlin(lib, fortran, i, xi, yi, nxi, nyi):
    psi, psi_ni = c_double(), c_double()
    if fortran:
        lib.__m_xpanel_MOD_pswlin(byref(c_int(i + 1)), byref(c_double(xi)), byref(c_double(yi)),
                                  byref(c_double(nxi)), byref(c_double(nyi)), byref(psi), byref(psi_ni))
    else:
        lib.pswlin(c_int(i), c_double(xi), c_double(yi), c_double(nxi), c_double(nyi), byref(psi), byref(psi_ni))
    return psi.value.hex(), psi_ni.value.hex()


def points(st, rng, n, nw, n_off=15):
    """(i, x, y, nx, ny): every airfoil node, every wake node, and off-body points (i = -1)."""
    xs, _ = st.views('x'); ys, _ = st.views('y'); nxs, _ = st.views('nx'); nys, _ = st.views('ny')
    pts = [(i, xs[i], ys[i], nxs[i], nys[i]) for i in range(n + nw)]
    for _ in range(n_off):
        a = rng.uniform(0, 2 * np.pi)
        pts.append((-1, float(rng.uniform(-1, 2)), float(rng.uniform(-1, 1)), math.cos(a), math.sin(a)))
    return pts


STATES = [(seed, closed, sharp) for seed in range(8) for closed, sharp in ((False, False), (True, False), (False, True))]


@pytest.mark.parametrize('seed,closed_te,sharp_le', STATES)
def test_apcalc_matches_fortran(clib, flib, seed, closed_te, sharp_le):
    st, _, _, _ = setup_state(clib, flib, seed, closed_te, sharp_le)
    clib.apcalc()
    flib.__m_xpanel_MOD_apcalc()
    assert st.same('apanel')


@pytest.mark.parametrize('seed,closed_te,sharp_le', STATES)
def test_psilin_matches_fortran(clib, flib, seed, closed_te, sharp_le):
    st, rng, n, nw = setup_state(clib, flib, seed, closed_te, sharp_le)
    clib.apcalc()
    flib.__m_xpanel_MOD_apcalc()
    for limage in (0, 1):
        st.set('limage', limage, np.int32)
        for k, (i, xi, yi, nxi, nyi) in enumerate(points(st, rng, n, nw)):
            combos = [(0, 0), (0, 1), (1, 0), (1, 1)] if (i < 0 or k % 7 == 0) else [(k % 2, (k // 2) % 2)]
            for geolin, siglin in combos:
                c = call_psilin(clib, False, i, xi, yi, nxi, nyi, geolin, siglin)
                f = call_psilin(flib, True, i, xi, yi, nxi, nyi, geolin, siglin)
                assert c == f, (limage, i, geolin, siglin)
                for name in WRITTEN:
                    assert st.same(name), (name, limage, i, geolin, siglin)


@pytest.mark.parametrize('seed,closed_te,sharp_le', STATES)
def test_pswlin_matches_fortran(clib, flib, seed, closed_te, sharp_le):
    st, rng, n, nw = setup_state(clib, flib, seed, closed_te, sharp_le)
    clib.apcalc()
    flib.__m_xpanel_MOD_apcalc()
    apc, apf = st.views('apanel')
    wake_ang = rng.uniform(-np.pi, np.pi, nw)          # wake panel angles (set by xywake in XFOIL)
    apc[n:n + nw] = wake_ang
    apf[n:n + nw] = wake_ang
    for i, xi, yi, nxi, nyi in points(st, rng, n, nw):
        assert call_pswlin(clib, False, i, xi, yi, nxi, nyi) == call_pswlin(flib, True, i, xi, yi, nxi, nyi), i
        for name in WRITTEN:
            assert st.same(name), (name, i)


def test_ncalc_matches_fortran(clib, flib):
    corners = 0
    for _, x, y, s, _, _ in G.shapes(clib, n=100):
        n = len(x)
        out = []
        for fortran in (False, True):
            xn, yn = np.zeros(n), np.zeros(n)
            if fortran:
                flib.__m_xpanel_MOD_ncalc(ptr(x), ptr(y), ptr(s), byref(c_int(n)), ptr(xn), ptr(yn))
            else:
                clib.ncalc(ptr(x), ptr(y), ptr(s), c_int(n), ptr(xn), ptr(yn))
            out.append(xn.tobytes() + yn.tobytes())
        assert out[0] == out[1]
        corners += int(np.any(np.diff(s) == 0))
    assert corners > 10


@pytest.mark.parametrize('seed,closed_te,sharp_le', STATES)
def test_ggcalc_qwcalc_qdcalc_match_fortran(clib, flib, seed, closed_te, sharp_le):
    """Whole-state comparison (all 428 globals) after each step of the panel solution."""
    from conftest import globals_diff, sync_globals
    sync_globals(clib, flib)
    st, rng, n, nw = setup_state(clib, flib, seed, closed_te, sharp_le)
    for flag in ('lqaij', 'lgamu', 'ladij', 'lwdij', 'limage'):
        st.set(flag, 0, np.int32)
    clib.apcalc()
    flib.__m_xpanel_MOD_apcalc()
    apc, apf = st.views('apanel')
    wake_ang = rng.uniform(-np.pi, np.pi, nw)
    apc[n:n + nw] = wake_ang
    apf[n:n + nw] = wake_ang
    assert globals_diff(clib, flib) == []
    for step in ('ggcalc', 'qwcalc', 'qdcalc', 'qdcalc'):        # 2nd qdcalc: ladij already set
        getattr(clib, step)()
        getattr(flib, '__m_xpanel_MOD_' + step)()
        assert globals_diff(clib, flib) == [], step
    c_piv, _ = st.views('aijpiv', np.int32)
    assert c_piv[:n + 1].min() >= 0 and c_piv[:n + 1].max() <= n         # 0-based pivot rows
    flags = {f: st.views(f, np.int32)[0][0] for f in ('lqaij', 'lgamu', 'ladij', 'lwdij')}
    assert flags == {'lqaij': 1, 'lgamu': 1, 'ladij': 1, 'lwdij': 1}


def run_wake_chain(clib, flib, seed, closed_te, sharp_le):
    """apcalc -> ggcalc -> gam(alfa) -> xywake -> stfind -> iblpan -> xicalc, whole state after each."""
    from conftest import globals_diff, sync_globals
    sync_globals(clib, flib)
    st, rng, n, _ = setup_state(clib, flib, seed, closed_te, sharp_le)
    for flag in ('lqaij', 'lgamu', 'ladij', 'lwdij', 'limage', 'lwake', 'lipan'):
        st.set(flag, 0, np.int32)
    waklen = [0.5, 1.0, 2.5, 6.0][seed % 4]                # 6.0: nw = n/12 + 60 > IWX for large n
    st.set('waklen', waklen)
    x, _ = st.views('x')
    st.set('chord', float(np.ptp(x[:n])))
    clib.apcalc()
    flib.__m_xpanel_MOD_apcalc()
    clib.ggcalc()
    flib.__m_xpanel_MOD_ggcalc()
    alfa = st.views('alfa')[0][0]
    gamu = st.views('gamu')[0].reshape((IQX, 2), order='F')
    st.set('gam', math.cos(alfa) * gamu[:n, 0] + math.sin(alfa) * gamu[:n, 1])
    assert globals_diff(clib, flib) == []
    for step in ('xywake', 'stfind', 'iblpan', 'xicalc'):
        getattr(clib, step)()
        getattr(flib, '__m_xpanel_MOD_' + step)()
        assert globals_diff(clib, flib) == [], step
    return st, n


@pytest.mark.parametrize('seed,closed_te,sharp_le', STATES)
def test_wake_and_bl_mapping_match_fortran(clib, flib, seed, closed_te, sharp_le):
    st, n = run_wake_chain(clib, flib, seed, closed_te, sharp_le)
    nw = st.views('nw', np.int32)[0][0]
    ist = st.views('ist', np.int32)[0][0]
    iblte = st.views('iblte', np.int32)[0][:2]
    nbl = st.views('nbl', np.int32)[0][:2]
    assert 0 <= ist < n - 1                                   # 0-based stagnation node
    assert iblte[0] == ist + 1 and iblte[1] == n - 1 - ist    # 0-based TE stations
    assert nbl[0] == iblte[0] + 1 and nbl[1] == iblte[1] + nw + 1
    ipan = st.views('ipan', np.int32)[0].reshape((452, 2), order='F')
    assert ipan[1, 0] == ist and ipan[iblte[0], 0] == 0 and ipan[nbl[1] - 1, 1] == n + nw - 1


def test_xywake_iwx_limit_reached(clib, flib):
    capped = 0
    for seed in range(3, 24, 4):                               # waklen = 6.0
        st, n = run_wake_chain(clib, flib, seed, False, False)
        capped += st.views('nw', np.int32)[0][0] == 82         # IWX
    assert capped > 0


def test_stfind_fallback_matches_fortran(clib, flib):
    from conftest import globals_diff, sync_globals
    sync_globals(clib, flib)
    st, rng, n, _ = setup_state(clib, flib, 0)
    st.set('gam', np.abs(rng.normal(size=n)) + 0.1)           # no sign change: "not found"
    clib.stfind()
    flib.__m_xpanel_MOD_stfind()
    assert globals_diff(clib, flib) == []
    assert st.views('ist', np.int32)[0][0] == n // 2 - 1


IBLPAN_STOP = '''
import sys
from ctypes import CDLL, byref, c_bool, c_char, c_size_t, c_void_p
import numpy as np
lib = CDLL(sys.argv[1]); fortran = sys.argv[2] == 'f'
if fortran:
    lib.set_print(byref(c_bool(True)))
else:
    lib.xfoil_default.restype = c_void_p
    lib.set_print(c_void_p(lib.xfoil_default()), byref(c_bool(True)))     # API: context first
def g(name):
    if fortran:
        return np.frombuffer((c_char * 4).in_dll(lib, '__i_xfoil_MOD_' + name), np.int32)
    lib.xfoil_global.restype = c_void_p
    return np.frombuffer((c_char * 4).from_address(lib.xfoil_global(name.encode(), None)), np.int32)
g('n')[0] = 600; g('nw')[0] = 82
g('ist')[0] = 590 if fortran else 589          # same stagnation node (1-based / 0-based)
if fortran:
    lib.__m_xpanel_MOD_iblpan()
else:
    lib.iblpan(c_void_p(lib.xfoil_default()))      # internal routine: context first
print('NOT STOPPED')
'''


def test_iblpan_overflow_stop_matches_fortran(clib, flib):
    import subprocess, sys
    from conftest import FORTRAN_LIB
    from xfoil_c import _lib
    res = {}
    for tag, path in (('c', _lib.lib_path), ('f', FORTRAN_LIB)):
        p = subprocess.run([sys.executable, '-c', IBLPAN_STOP, path, tag], capture_output=True, text=True)
        res[tag] = (p.returncode, p.stdout, p.stderr)
    assert res['c'][:2] == res['f'][:2]                 # exit status 0 and the same printed lines
    assert res['c'][0] == 0 and 'Increase IVX to at least         673' in res['c'][1]
    # Fortran has a bare `stop` (nothing on stderr); C gives it a message since clean-up C2, so that
    # Python's XFoilError can say why (outside an API call it is printed like gfortran's STOP text)
    assert res['f'][2] == '' and res['c'][2] == 'STOP IBLPAN: BL array overflow (increase IVX)\n'


IVX = 452
BL_ARRAYS = ('ctau', 'thet', 'dstr', 'uedg', 'mass')


@pytest.mark.parametrize('seed,closed_te,sharp_le', STATES[::2])
def test_viscous_coupling_routines_match_fortran(clib, flib, seed, closed_te, sharp_le):
    """qwcalc, qdcalc, qiset, uicalc, qvfue, gamqv, ueset, stmove (stagnation point moved
    up, down and not at all); all 428 globals identical after every call."""
    from conftest import globals_diff
    st, n = run_wake_chain(clib, flib, seed, closed_te, sharp_le)

    def step(name):
        getattr(clib, name)()
        getattr(flib, '__m_xpanel_MOD_' + name)()
        assert globals_diff(clib, flib) == [], name

    for name in ('qwcalc', 'qdcalc', 'qiset', 'uicalc'):
        step(name)

    rng = np.random.default_rng(1000 + seed)
    for arr in BL_ARRAYS:
        st.set(arr, rng.uniform(0.01, 1.0, IVX * 2))
    iblte_c = st.views('iblte', np.int32)[0][:2]
    st.set_index('itran', [int(rng.integers(2, iblte_c[0] + 2)), int(rng.integers(2, iblte_c[1] + 2))])
    assert globals_diff(clib, flib) == []

    for name in ('qvfue', 'gamqv', 'ueset'):
        step(name)

    moved = []
    for shift in (0, 3, -4, 1, -1):
        ist_old = int(st.views('ist', np.int32)[0][0])
        m = min(max(ist_old + shift, 2), n - 4)
        st.set('gam', (m + 0.5 - np.arange(n)) * 0.01)        # sign change between nodes m, m+1
        uc, uf = st.views('uedg')
        k = int(rng.integers(1, 5))
        uc[k], uf[k] = -1e-3, -1e-3                              # below ueps: exercises the tweak
        step('stmove')
        ist_new = int(st.views('ist', np.int32)[0][0])
        assert ist_new == m
        moved.append(ist_new - ist_old)
    assert moved[0] == 0 and any(d > 0 for d in moved) and any(d < 0 for d in moved)
