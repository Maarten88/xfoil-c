"""Phase 10b flap + hinge moment: C port of the upstream XFOIL 6.99 FLAP vs that routine itself
(tests/flapref: verbatim upstream FLAP/INSIDE/GETXYF/SSS/SCHECK, built here with the reference
flags and linked to libxfoil_fortran.so for the spline/geometry routines), then end to end
against the Fortran library (hinge moments through mhinge), and the Python API."""
import shutil
import subprocess
from ctypes import CDLL, byref, c_bool, c_double, c_int, c_void_p
from pathlib import Path

import numpy as np
import pytest

import xfoil_c
from conftest import FORTRAN_LIB, global_array, globals_diff, sync_globals
import test_xfoil as XF
from test_api_analysis import naca0012

HERE = Path(__file__).resolve().parent / 'flapref'
FFLAGS = ['-O2', '-fdefault-real-8', '-fno-range-check', '-fno-tree-vectorize', '-std=legacy']


@pytest.fixture(scope='session')
def ref(tmp_path_factory, flib):
    if not shutil.which('gfortran'):
        pytest.skip('gfortran not available')
    d = tmp_path_factory.mktemp('flapref')
    libdir = str(Path(FORTRAN_LIB).parent)
    subprocess.run(['gcc', '-c', '-fPIC', '-O2', '-ffp-contract=off', str(HERE / 'glue.c'), '-o', str(d / 'glue.o')],
                   check=True)
    subprocess.run(['gfortran', '-shared', '-fPIC', *FFLAGS, '-I' + str(HERE),
                    *[str(HERE / f) for f in ('upstream_xgdes.f', 'upstream_xgeom.f', 'driver.f')],
                    str(d / 'glue.o'), '-o', str(d / 'libflapref.so'),
                    '-L' + libdir, '-l:' + Path(FORTRAN_LIB).name, '-Wl,-rpath,' + libdir],
                   check=True, capture_output=True)
    return CDLL(str(d / 'libflapref.so'))


def ptr(a):
    return a.ctypes.data_as(c_void_p)


def ref_flap(ref, x, y, xh, yh, ddef, lgsym=False):
    """Upstream FLAP on the buffer x, y (splined by the driver like set_airfoil does)."""
    n = len(x)
    out = {k: np.zeros(4 * 370) for k in ('xb', 'yb', 'sb', 'xbp', 'ybp')}
    pars = np.zeros(12)
    nout = c_int()
    ref.refflap_(byref(c_int(n)), ptr(np.array(x, float)), ptr(np.array(y, float)), byref(c_double(xh)),
                 byref(c_double(yh)), byref(c_double(ddef)), byref(c_int(int(lgsym))), byref(nout),
                 *[ptr(out[k]) for k in ('xb', 'yb', 'sb', 'xbp', 'ybp')], ptr(pars))
    m = nout.value
    res = {k: v[:m].copy() for k, v in out.items()}
    res.update(nb=m, **dict(zip(('sble', 'chordb', 'areab', 'radble', 'angbte', 'thickb', 'cambrb', 'xbf', 'ybf',
                                 'lbflap', 'lgsame', 'lgsym'), pars)))
    return res


def c_buffer(clib):
    nb = int(global_array(clib, 'nb', np.int32)[0])
    size = 4 * 640
    res = {k: np.array(global_array(clib, k, np.float64, (size,), False)[:nb]) for k in ('xb', 'yb', 'sb', 'xbp', 'ybp')}
    res['nb'] = nb
    for k in ('sble', 'chordb', 'areab', 'radble', 'angbte', 'thickb', 'cambrb', 'xbf', 'ybf'):
        res[k] = float(global_array(clib, k, np.float64)[0])
    for k in ('lbflap', 'lgsym'):                     # (lgsame: dropped in C, clean-up C3)
        res[k] = float(global_array(clib, k, np.int32)[0])
    return res


def same(a, b):
    assert a['nb'] == b['nb'], (a['nb'], b['nb'])
    for k in a:
        assert np.array(a[k], float).tobytes() == np.array(b[k], float).tobytes(), k


def hinge_y(lib, fortran, x_hinge, y_rel, flib):
    """y of the hinge for y/t = y_rel at x_hinge on lib's buffer airfoil, as GETXYF computes it
    (Fortran spline routines of flib)."""
    nb = int(global_array(lib, 'nb', np.int32, None, fortran)[0])
    arr = {k: np.array(global_array(lib, k, np.float64, (4 * 640,), fortran)) for k in ('xb', 'yb', 'sb', 'xbp', 'ybp')}
    seval = flib_fn(flib, 'seval')
    sinvrt = flib_fn(flib, 'sinvrt')
    tops = c_double(arr['sb'][0] + (arr['xb'][0] - x_hinge))
    bots = c_double(arr['sb'][nb - 1] - (arr['xb'][nb - 1] - x_hinge))
    for sv in (tops, bots):
        sinvrt(byref(sv), byref(c_double(x_hinge)), ptr(arr['xb']), ptr(arr['xbp']), ptr(arr['sb']), byref(c_int(nb)))
    topy = seval(byref(tops), ptr(arr['yb']), ptr(arr['ybp']), ptr(arr['sb']), byref(c_int(nb)))
    boty = seval(byref(bots), ptr(arr['yb']), ptr(arr['ybp']), ptr(arr['sb']), byref(c_int(nb)))
    return topy * y_rel + boty * (1.0 - y_rel)


def flib_fn(flib, name):
    f = getattr(flib, '__m_spline_MOD_' + name)
    f.restype = c_double if name == 'seval' else None
    return f


def start(clib, flib, seed):
    """Both libraries: init + the same airfoil (NACA or random shape), maybe repaneled."""
    rng = np.random.default_rng(700 + seed)
    XF.init_both(clib, flib, rng)
    if seed % 3 == 0:
        x, y = naca0012()
    elif seed % 3 == 1:
        x, y = XF.buffer_shape(rng, False, 0)
    else:
        spec = [2412, 4415, 23012][(seed // 3) % 3]
        for lib in (clib, flib):
            lib.set_naca(byref(c_int(spec)))
        x = None
    if x is not None:
        for lib in (clib, flib):
            lib.set_airfoil(ptr(np.array(x, float)), ptr(np.array(y, float)), byref(c_int(len(x))))
    assert globals_diff(clib, flib) == []
    return rng


CASES = [(s, a) for s in range(8) for a in (5.0, -7.5, 15.0, 0.4, -25.0)]


@pytest.mark.parametrize('seed,angle', CASES)
def test_flap_buffer_matches_upstream(clib, flib, ref, seed, angle):
    rng = start(clib, flib, seed)
    x_hinge = float(rng.uniform(0.55, 0.85))
    y_rel = float(rng.choice([0.5, rng.uniform(0.1, 0.9), -0.3, 1.4]))    # inside, or below/above the airfoil
    before = c_buffer(clib)
    clib.set_hinge(byref(c_double(x_hinge)), byref(c_double(y_rel)), byref(c_bool(True)))
    yh = float(global_array(clib, 'ybf', np.float64)[0])
    assert yh == hinge_y(clib, False, x_hinge, y_rel, flib)
    clib.flap(byref(c_double(angle)))
    same(c_buffer(clib), ref_flap(ref, before['xb'], before['yb'], x_hinge, yh, angle))
    # absolute deflection (user decision 2026-11-21): the same angle again changes nothing, another
    # angle starts from the undeflected airfoil (XFOIL's FLAP would add up)
    once = c_buffer(clib)
    clib.flap(byref(c_double(angle)))
    same(c_buffer(clib), once)
    clib.flap(byref(c_double(-angle / 2)))
    same(c_buffer(clib), ref_flap(ref, before['xb'], before['yb'], x_hinge, yh, -angle / 2))
    clib.flap(byref(c_double(0.0)))
    restored = c_buffer(clib)
    assert restored['nb'] == before['nb']
    for k in ('xb', 'yb'):
        assert np.array(restored[k]).tobytes() == np.array(before[k]).tobytes(), k


def test_flap_coverage(clib, flib, ref):
    """The cases reach both break types (gap opened with new arc points / corner points) on both
    surfaces and add and remove points."""
    seen = set()
    for seed, angle in CASES:
        rng = start(clib, flib, seed)
        x_hinge = float(rng.uniform(0.55, 0.85))
        y_rel = float(rng.choice([0.5, rng.uniform(0.1, 0.9), -0.3, 1.4]))
        nb0 = int(global_array(clib, 'nb', np.int32)[0])
        clib.set_hinge(byref(c_double(x_hinge)), byref(c_double(y_rel)), byref(c_bool(True)))
        clib.flap(byref(c_double(angle)))
        nb1 = int(global_array(clib, 'nb', np.int32)[0])
        seen.add(('inside' if 0 <= y_rel <= 1 else 'outside', 'down' if angle > 0 else 'up'))
        seen.add('more points' if nb1 > nb0 else 'fewer points' if nb1 < nb0 else 'same count')
    assert {('inside', 'down'), ('inside', 'up'), ('outside', 'down'), ('outside', 'up'),
            'more points'} <= seen, seen


@pytest.mark.parametrize('seed', range(6))
def test_hinge_moment_end_to_end_matches_fortran(clib, flib, ref, seed):
    """C set_hinge/flap + analyses vs the Fortran library given the same hinge and the upstream
    deflected buffer (then Fortran pangen, as flap does), whole state and outputs compared."""
    rng = start(clib, flib, seed)
    x_hinge = float(rng.uniform(0.6, 0.8))
    y_rel = float(rng.uniform(0.2, 0.8))
    angle = [0.0, 6.0, -4.0][seed % 3]
    for lib in (clib, flib):
        lib.set_reynolds(byref(c_double([0.0, 5e5, 1e6][seed % 3])))
        lib.set_max_iter(byref(c_int(60)))
    before = c_buffer(clib)
    clib.set_hinge(byref(c_double(x_hinge)), byref(c_double(y_rel)), byref(c_bool(True)))
    yh = hinge_y(flib, True, x_hinge, y_rel, flib)
    # Fortran side: the same hinge (FMOM/GETXYF path) ...
    for name, v in (('xbf', x_hinge), ('ybf', yh), ('xof', x_hinge), ('yof', yh),
                    ('hmom', np.nan), ('hfx', np.nan), ('hfy', np.nan)):
        global_array(flib, name, np.float64, None, True)[:] = v
    for name in ('lbflap', 'lflap'):
        global_array(flib, name, np.int32, None, True)[:] = 1
    assert globals_diff(clib, flib) == []
    if angle != 0.0:
        clib.flap(byref(c_double(angle)))
        # ... and the upstream deflection, then geopar + pangen as the C flap does
        r = ref_flap(ref, before['xb'], before['yb'], x_hinge, yh, angle)
        nb = r['nb']
        for k in ('xb', 'yb', 'sb', 'xbp', 'ybp'):
            fa = global_array(flib, k, np.float64, (4 * 640,), True)
            fa[:nb] = r[k]
            # beyond nb: stale values left by FLAP's array shifts (unused); take them from C
            fa[nb:] = global_array(clib, k, np.float64, (4 * 640,), False)[nb:]
            assert np.array(global_array(clib, k, np.float64, (4 * 640,), False)[:nb]).tobytes() == r[k].tobytes()
        global_array(flib, 'nb', np.int32, None, True)[:] = nb
        global_array(flib, 'lgsame', np.int32, None, True)[:] = 0
        a = {k: global_array(flib, k, np.float64, (4 * 640,), True) for k in ('xb', 'xbp', 'yb', 'ybp', 'sb', 'w1')}
        outs = {k: global_array(flib, k, np.float64, None, True) for k in
                ('sble', 'chordb', 'areab', 'radble', 'angbte', 'ei11ba', 'ei22ba', 'apx1ba', 'apx2ba',
                 'ei11bt', 'ei22bt', 'apx1bt', 'apx2bt', 'thickb', 'cambrb')}
        flib.__m_xgeom_MOD_geopar(*[ptr(a[k]) for k in ('xb', 'xbp', 'yb', 'ybp', 'sb')], byref(c_int(nb)),
                                  ptr(a['w1']), *[ptr(v) for v in outs.values()])
        flib.__m_xfoil_MOD_pangen(byref(c_int(1)))
        for name in ('hmom', 'hfx', 'hfy'):
            global_array(flib, name, np.float64, None, True)[:] = np.nan
        assert globals_diff(clib, flib) == []
    for alpha in (2.0, 6.0):
        res = []
        for lib in (clib, flib):
            out = [c_double() for _ in range(6)]
            conv = c_bool()
            lib.alfa(byref(c_double(alpha)), *[byref(o) for o in out], byref(conv))
            res.append([o.value for o in out] + [conv.value])
        assert np.array(res[0]).tobytes() == np.array(res[1]).tobytes()
        assert globals_diff(clib, flib) == []
        h = [c_double() for _ in range(3)]
        clib.get_hinge(*[byref(v) for v in h])
        assert [v.value for v in h] == [float(global_array(flib, k, np.float64, None, True)[0])
                                        for k in ('hmom', 'hfx', 'hfy')]
        assert np.isfinite(h[0].value)


def test_python_api(clib):
    xf = xfoil_c.XFoil()
    xf.print = False
    global_array(clib, 'nb', np.int32)[:] = 0              # as in a fresh process (init keeps the buffer)
    assert np.isnan(xf.hinge_moment) and xf.hinge is None
    with pytest.raises(ValueError):
        xf.set_hinge(0.75, 0.5, y_relative=True)            # no airfoil yet
    xf.naca('2412')
    with pytest.raises(RuntimeError):
        xf.flap(5.0)                                        # no hinge yet
    xf.Re = 1e6
    xf.max_iter = 100
    xf.set_hinge(0.75, 0.5, y_relative=True)
    assert xf.hinge[0] == 0.75 and np.isnan(xf.hinge_moment)
    cl0 = xf.a(4.0)[0]
    hm0 = xf.hinge_moment
    assert np.isfinite(hm0) and np.isfinite(xf.hinge_fx) and np.isfinite(xf.hinge_fy)
    n0 = xf.airfoil.n_coords
    x0 = xf.airfoil.x.copy()
    assert xf.flap_angle == 0.0
    xf.flap(0.0)                                            # undeflected (repanels with the same settings)
    assert xf.airfoil.n_coords == n0 and np.array_equal(xf.airfoil.x, x0)
    xf.flap(5.0)                                            # TE down: more lift, different hinge moment
    assert np.isnan(xf.hinge_moment) and xf.flap_angle == 5.0
    cl5 = xf.a(4.0)[0]
    assert cl5 > cl0 + 0.2 and abs(xf.hinge_moment - hm0) > 1e-3
    x5 = xf.airfoil.x.copy()
    xf.flap(5.0)                                            # absolute: the same angle again changes nothing
    assert np.array_equal(xf.airfoil.x, x5) and abs(xf.a(4.0)[0] - cl5) < 1e-3
    xf.flap(-5.0)                                           # TE up from the undeflected airfoil
    assert xf.a(4.0)[0] < cl0 - 0.2
    xf.flap(0.0)                                            # back to the original section
    assert np.array_equal(xf.airfoil.x, x0) and abs(xf.a(4.0)[0] - cl0) < 1e-3
    xf.flap(5.0)
    xf.set_hinge(0.65, 0.5, y_relative=True)                # moving the hinge keeps the deflection
    assert xf.flap_angle == 5.0 and xf.hinge[0] == 0.65 and not np.array_equal(xf.airfoil.x, x5)
    xf.airfoil = xf.airfoil                                 # a new airfoil clears the hinge and deflection
    assert xf.hinge is None and np.isnan(xf.hinge_moment) and xf.flap_angle == 0.0
    xf.set_hinge(0.7, 0.0)
    xf.naca('0012')                                         # so does naca()
    assert xf.hinge is None
    for name in ('hinge_moment', 'hinge_fx', 'hinge_fy', 'hinge', 'flap_angle'):
        with pytest.raises(AttributeError):
            setattr(xf, name, 1.0)
