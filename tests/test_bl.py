"""Extension get_n_bl/get_bl (boundary layer of the last analysed point): the C values against the file
the Fortran reference writes with its OPER DUMP routine (m_xoper bldump) after the same API calls (whole
state compared after each call); NaN rules; the XFoil.get_bl/get_bl_nodes wrappers."""
from ctypes import byref, c_bool, c_char_p, c_double, c_int, c_size_t, c_void_p

from pathlib import Path

import numpy as np
import pytest

from conftest import global_array, globals_diff, sync_globals
from test_api import call

FIELDS = ('s', 'x', 'y', 'ue', 'dstar', 'theta', 'cf', 'hk', 'hstar', 'n', 'ctau')
# bldump format 99003: (1x, 4F9.5, 3F10.6, f10.4, f10.4, ...): field widths and decimals of the
# columns s x y ue dstar theta cf hk hstar; wake lines stop after hk
WIDTHS = (9, 9, 9, 9, 10, 10, 10, 10, 10)
DECIMALS = (5, 5, 5, 5, 6, 6, 6, 4, 4)

CASES = [  # naca, n_panel, Re, Mach, Ncrit, alpha
    ('0012', 160, 1e6, 0.0, 9.0, 4.0),
    ('4415', 160, 5e5, 0.0, 2.0, 8.0),
    ('2412', 200, 3e5, 0.3, 9.0, -2.0),      # compressible Hk (hstinv != 0)
    ('23012', 120, 2e6, 0.0, 5.0, 10.0),
]


def fortran_dump(flib, path):
    """path: a short relative name (bldump copies it into the fixed-length FNAME)."""
    c_int.in_dll(flib, '__i_xfoil_MOD_kdelim').value = 0            # blank-delimited, all columns
    name = str(path).encode()
    flib.__m_xoper_MOD_bldump(c_char_p(name), c_size_t(len(name)))
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            continue
        vals, pos = [], 1
        for w in WIDTHS:
            field = line[pos:pos + w]
            pos += w
            vals.append(float(field) if field.strip() else np.nan)
        rows.append(vals)
    return np.array(rows)


def c_bl(clib):
    n = clib.get_n_bl()
    out = [np.zeros(n) for _ in FIELDS]
    clib.get_bl(*[a.ctypes.data_as(c_void_p) for a in out], byref(c_int(n)))
    return dict(zip(FIELDS, out))


def analyse(clib, flib, naca, n_panel, re, mach, n_crit, alpha):
    sync_globals(clib, flib)
    clib.init(); flib.init()
    call(clib, flib, 'set_naca', int(naca))
    call(clib, flib, 'repanel', n_panel, 1.0, 0.15, 0.2, 1.0, 1.0, 1.0, 1.0)
    call(clib, flib, 'set_reynolds', re)
    call(clib, flib, 'set_mach', mach)
    call(clib, flib, 'set_n_crit', n_crit)
    call(clib, flib, 'set_max_iter', 200)
    conv = []
    for lib in (clib, flib):
        outs = [c_double() for _ in range(6)]
        ok = c_bool()
        lib.alfa(byref(c_double(alpha)), *[byref(o) for o in outs], byref(ok))
        conv.append(ok.value)
    assert globals_diff(clib, flib) == []
    return conv[0]


@pytest.mark.parametrize('case', CASES, ids=[f'naca{c[0]}_a{c[5]:g}' for c in CASES])
def test_get_bl_matches_fortran_bldump(clib, flib, tmp_path, monkeypatch, case):
    assert analyse(clib, flib, *case)
    bl = c_bl(clib)
    monkeypatch.chdir(tmp_path)
    ref = fortran_dump(flib, Path('bl.dat'))
    n = clib.get_n_cp()
    assert ref.shape[0] == len(bl['x']) == n + int(global_array(clib, 'nw', np.int32)[0])
    for k, (name, dec) in enumerate(zip(FIELDS[:9], DECIMALS)):
        c = bl[name].copy()
        if name == 'hstar':
            c[n:] = np.nan                         # bldump writes no H* in the wake
        r = ref[:, k]
        tol = 0.5 * 10.0 ** -dec + 1e-12
        both = np.isfinite(r) & np.isfinite(c)
        assert np.array_equal(np.isfinite(r), np.isfinite(c)), name
        assert np.max(np.abs(c[both] - r[both])) <= tol, name
    # s, x, y are the current-airfoil (and wake) nodes, bit for bit
    for name in ('s', 'x', 'y'):
        assert np.array_equal(bl[name], global_array(clib, name, np.float64)[:len(bl[name])]), name


def test_laminar_turbulent_split(clib, flib):
    assert analyse(clib, flib, *CASES[1])
    bl = c_bl(clib)
    lam, turb = np.isfinite(bl['n']), np.isfinite(bl['ctau'])
    assert not np.any(lam & turb) and np.all(lam | turb)
    n = clib.get_n_cp()
    assert np.all(turb[n:])                        # the wake is turbulent
    assert np.all(bl['cf'][n:] == 0.0)
    assert np.all(bl['n'][lam] >= 0.0) and np.all(bl['n'][lam] <= 2.0 + 1e-9)   # Ncrit 2


def test_nan_without_converged_viscous_point(clib, flib):
    sync_globals(clib, flib)
    clib.init()
    clib.set_naca(byref(c_int(12)))
    bl = c_bl(clib)                               # nothing analysed yet
    assert np.all(np.isfinite(bl['x'])) and np.all(np.isnan(bl['ue'])) and np.all(np.isnan(bl['cf']))
    clib.set_reynolds(byref(c_double(0.0)))       # inviscid point: converged, but no boundary layer
    outs = [c_double() for _ in range(6)]
    ok = c_bool()
    clib.alfa(byref(c_double(3.0)), *[byref(o) for o in outs], byref(ok))
    assert ok.value
    assert np.all(np.isnan(c_bl(clib)['theta']))


def test_wrappers():
    from xfoil_c import XFoil
    xf = XFoil()
    xf.print = False
    xf.naca('4415')
    xf.repanel()
    xf.Re = 5e5
    xf.n_crit = 2.0
    xf.max_iter = 200
    assert xf.get_bl() is None
    cl, cd, cm, cp, xtr_top, xtr_bot = xf.a(6.0)
    nodes = xf.get_bl_nodes()
    assert set(nodes) == set(XFoil.BL_FIELDS)
    bl = xf.get_bl()
    n_air = len(xf.get_cp_distribution()[0])
    assert len(bl['top']['x']) + len(bl['bottom']['x']) == n_air
    assert len(bl['wake']['x']) == len(nodes['x']) - n_air
    for side in ('top', 'bottom'):
        d = bl[side]
        assert np.all(d['ue'] >= 0.0)
        assert d['x'][-1] == pytest.approx(1.0, abs=1e-3)             # ends at the trailing edge
        turb = np.nonzero(np.isfinite(d['ctau']))[0]
        assert np.all(np.isfinite(d['n'][:turb[0]]))                  # laminar up to transition
    # the first turbulent top station lies just downstream of XFoil's transition point
    x_first = bl['top']['x'][np.nonzero(np.isfinite(bl['top']['ctau']))[0][0]]
    assert xtr_top <= x_first < xtr_top + 0.05
    xf.a(40.0)                                     # does not converge
    assert xf.get_bl() is None
