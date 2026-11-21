"""Phase 10b extensions: get_geometry (C) vs the Fortran geopar/tccalc on the same current airfoil;
it must not change the state; the XFoil properties report the same numbers."""
from ctypes import byref, c_double, c_int, c_void_p

import numpy as np
import pytest

import xfoil_c
from conftest import global_array, globals_diff, sync_globals
import test_xfoil as XF


def _izx():
    """Length of the current-airfoil arrays (x, y, ...: IZX = IQX + IWX) from the generated globals."""
    from ctypes import c_size_t
    from xfoil_c import _lib
    size = c_size_t()
    _lib.lib.xfoil_global(b'x', byref(size))
    return size.value // 8


IZX = _izx()


def current_airfoil(lib, fortran):
    n = int(global_array(lib, 'n', np.int32, None, fortran)[0])
    arrs = {k: np.array(global_array(lib, k, np.float64, (IZX,), fortran)) for k in ('x', 'xp', 'y', 'yp', 's')}
    return n, arrs


def fortran_geometry(flib):
    n, a = current_airfoil(flib, True)
    t = np.zeros(IZX)
    outs = [c_double() for _ in range(15)]        # sle chord area radle angte ei11a .. apx2t thick cambr
    p = lambda k: a[k].ctypes.data_as(c_void_p)
    flib.__m_xgeom_MOD_geopar(p('x'), p('xp'), p('y'), p('yp'), p('s'), byref(c_int(n)),
                              t.ctypes.data_as(c_void_p), *[byref(o) for o in outs])
    thick, xthick, cambr, xcambr = (c_double() for _ in range(4))
    flib.__m_xgeom_MOD_tccalc(p('x'), p('xp'), p('y'), p('yp'), p('s'), byref(c_int(n)),
                              byref(thick), byref(xthick), byref(cambr), byref(xcambr))
    o = [v.value for v in outs]
    return [o[13], xthick.value, o[14], xcambr.value, o[3], o[4], o[2]]


def c_geometry(clib):
    out = [c_double() for _ in range(7)]
    clib.get_geometry(*[byref(v) for v in out])
    return [v.value for v in out]


@pytest.mark.parametrize('seed', range(10))
def test_get_geometry_matches_fortran_and_keeps_state(clib, flib, seed):
    rng = np.random.default_rng(500 + seed)
    XF.init_both(clib, flib, rng)
    if seed < 4:
        spec = [12, 2412, 4415, 23012][seed]
        clib.set_naca(byref(c_int(spec))); flib.set_naca(byref(c_int(spec)))
    else:
        x, y = XF.buffer_shape(rng, bool(seed % 2), seed % 3)
        for lib in (clib, flib):
            lib.set_airfoil(x.ctypes.data_as(c_void_p), y.ctypes.data_as(c_void_p), byref(c_int(len(x))))
        if seed % 2 == 0:
            npan = int(rng.integers(80, 250))          # drawn once: same panel count for both libraries
            for lib in (clib, flib):
                lib.repanel(byref(c_int(npan)), *[byref(c_double(v)) for v in (1.0, 0.15, 0.2, 1, 1, 1, 1)])
    assert globals_diff(clib, flib) == []
    c = c_geometry(clib)
    f = fortran_geometry(flib)
    assert np.array(c).tobytes() == np.array(f).tobytes(), (c, f)
    assert globals_diff(clib, flib) == []           # get_geometry did not change the C state
    assert c[0] > 0 and c[6] > 0


def test_get_geometry_without_airfoil(clib, flib):
    sync_globals(clib, flib)
    clib.init()
    XF.set_both(clib, flib, 'n', 0, np.int32)
    assert np.isnan(c_geometry(clib)).all()


def test_xfoil_properties(clib):
    xf = xfoil_c.XFoil()
    xf.print = False
    assert np.isnan(xf.max_thickness)               # init: no airfoil
    xf.naca('2412')
    g = c_geometry(xf._lib)                         # the instance's own context (C5b)
    assert (xf.max_thickness, xf.x_max_thickness, xf.max_camber, xf.x_max_camber, xf.le_radius, xf.area) == \
           (g[0], g[1], g[2], g[3], g[4], g[6])
    assert xf.te_angle == float(np.degrees(g[5]))
    assert abs(xf.max_thickness - 0.12) < 1e-3 and abs(xf.max_camber - 0.02) < 1e-3
    assert abs(xf.x_max_camber - 0.4) < 0.01
    for name in ('max_thickness', 'x_max_thickness', 'max_camber', 'x_max_camber', 'le_radius', 'te_angle', 'area'):
        with pytest.raises(AttributeError):
            setattr(xf, name, 1.0)                  # read-only
