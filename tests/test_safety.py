# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""Clean-up C6a: inputs that made the Fortran code read or write outside its arrays now behave
safely (and the library stays usable afterwards)."""
from ctypes import byref, c_double, c_int, c_void_p

import numpy as np
import pytest

import xfoil_c
from xfoil_c.model import Airfoil


def fresh():
    xf = xfoil_c.XFoil()
    xf.print = False
    return xf


def cang(clib, x, y, pad):
    """cang on x, y placed inside larger arrays whose extra elements hold `pad`."""
    n = len(x)
    xp = np.full(n + 4, pad); yp = np.full(n + 4, pad)
    xp[2:n + 2] = x; yp[2:n + 2] = y
    imax, amax = c_int(), c_double()
    clib.cang(xp[2:].ctypes.data_as(c_void_p), yp[2:].ctypes.data_as(c_void_p), c_int(n), c_int(0),
              byref(imax), byref(amax))
    return imax.value, amax.value


def test_cang_doubled_end_points_stay_inside_the_arrays(clib):
    t = np.linspace(0.0, 2.0 * np.pi, 41)
    x, y = 0.5 + 0.5 * np.cos(t), 0.1 * np.sin(t)
    x[1], y[1] = x[0], y[0]                      # doubled first point
    x[-2], y[-2] = x[-1], y[-1]                  # doubled last point
    results = {cang(clib, x, y, pad) for pad in (-7.0, 0.0, 3.0, np.nan)}
    assert len(results) == 1                     # independent of anything outside x[0..n-1]
    imax, amax = results.pop()
    assert 0 < imax < len(x) - 1 and np.isfinite(amax)


def test_flap_with_hinge_outside_the_airfoil_raises():
    xf = fresh()
    xf.naca('0012')
    xf.set_hinge(1.5, 0.0)                       # behind the TE: no surface breaks
    with pytest.raises(xfoil_c.XFoilError, match='FLAP: surface break points not found'):
        xf.flap(10.0)
    xf.naca('0012')                              # usable afterwards
    assert [round(v, 4) for v in xf.a(0.0)[:1]] == [0.0]


@pytest.mark.parametrize('n,match', [(3000, 'too many coordinate points'), (2, 'at least 3')])
def test_set_airfoil_point_count_is_checked(n, match):
    xf = fresh()
    t = np.linspace(0.0, 2.0 * np.pi, n)
    with pytest.raises(xfoil_c.XFoilError, match=match):
        xf.airfoil = Airfoil(0.5 + 0.5 * np.cos(t), 0.06 * np.sin(t))
    xf.naca('0012')
    assert np.isfinite(xf.a(2.0)[0])
