# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""m_naca.f90 -> src/xnaca.c: coordinates and names must be identical to the Fortran routines."""
from ctypes import byref, c_int, c_size_t, create_string_buffer, c_void_p

import numpy as np
import pytest

NSIDE = 640 // 3        # naca() uses IQX/3 points per side (IQX = 640)
NAME_LEN = 48           # i_xfoil: character(48) :: name


def ptr(arr):
    return arr.ctypes.data_as(c_void_p)   # typed pointer: a plain int is passed as 32-bit C int


def run(lib, routine, ides, nside, fortran):
    xx, yt, yc = np.zeros(nside), np.zeros(nside), np.zeros(nside)
    xb, yb = np.full(2 * nside, -1.0), np.full(2 * nside, -1.0)
    nb = c_int(-1)
    name = create_string_buffer(b'#' * NAME_LEN, NAME_LEN)
    if fortran:
        fn = getattr(lib, '__m_naca_MOD_' + routine)
        di = c_int(ides)
        fn(byref(di), ptr(xx), ptr(yt), ptr(yc), byref(c_int(nside)), ptr(xb), ptr(yb), byref(nb),
           name, c_size_t(NAME_LEN))
    else:
        if routine == 'naca4':
            di = c_int(ides)
            lib.naca4(di, ptr(xx), ptr(yt), ptr(yc), c_int(nside), ptr(xb), ptr(yb), byref(nb),
                      name, c_int(NAME_LEN))
        else:
            di = c_int(ides)
            lib.naca5(byref(di), ptr(xx), ptr(yt), ptr(yc), c_int(nside), ptr(xb), ptr(yb), byref(nb),
                      name, c_int(NAME_LEN))
    return di.value, nb.value, name.raw, xx, yt, yc, xb, yb


def assert_same(c, f):
    assert c[:3] == f[:3]                                   # ides, nb, name bytes
    for a, b in zip(c[3:], f[3:]):
        assert a.tobytes() == b.tobytes()


def test_naca4_all_designations(clib, flib):
    for ides in range(10000):
        assert_same(run(clib, 'naca4', ides, NSIDE, False), run(flib, 'naca4', ides, NSIDE, True))


@pytest.mark.parametrize('nside', [2, 3, 10, 57, 100])
def test_naca4_other_sizes(clib, flib, nside):
    for ides in (12, 2412, 4415, 9999, 1000, 10):
        assert_same(run(clib, 'naca4', ides, nside, False), run(flib, 'naca4', ides, nside, True))


def test_naca5_all_legal_designations(clib, flib):
    for n543 in (210, 220, 230, 240, 250):
        for tt in range(100):
            ides = n543 * 100 + tt
            c = run(clib, 'naca5', ides, NSIDE, False)
            assert_same(c, run(flib, 'naca5', ides, NSIDE, True))
            assert c[0] == ides and c[1] == 2 * NSIDE - 1


@pytest.mark.parametrize('ides', [0, 10012, 20012, 21112, 23512, 26012, 25599, 99999, 12345])
def test_naca5_illegal_designations(clib, flib, ides):
    c = run(clib, 'naca5', ides, NSIDE, False)
    assert_same(c, run(flib, 'naca5', ides, NSIDE, True))
    assert c[0] == 0 and c[1] == -1 and c[2] == b'#' * NAME_LEN     # outputs untouched


def test_naca_names(clib, flib):
    assert run(clib, 'naca4', 12, NSIDE, False)[2] == b'NACA 0012'.ljust(NAME_LEN)
    assert run(clib, 'naca5', 23012, NSIDE, False)[2] == b'NACA 23012'.ljust(NAME_LEN)


def test_naca4_0012_geometry(clib):
    _, nb, _, _, _, _, xb, yb = run(clib, 'naca4', 12, NSIDE, False)
    x, y = xb[:nb], yb[:nb]
    assert (x[0], x[-1], x[NSIDE - 1]) == (1.0, 1.0, 0.0)         # TE, TE, LE
    assert np.allclose(y, -y[::-1])                                 # symmetric
    assert abs(y.max() - 0.06) < 1e-3                               # 12% thick
