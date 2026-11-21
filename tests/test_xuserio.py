# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""m_userio.f90 strip -> src/xuserio.c: results must be byte-identical to the Fortran routine.

Fortran CHARACTER(*) arguments get a hidden trailing length argument (size_t in gfortran >= 8).
"""
from ctypes import byref, c_int, c_size_t, create_string_buffer

import numpy as np
import pytest


def run_strip(clib, flib, text):
    raw = text.encode('ascii')
    cbuf, fbuf = create_string_buffer(raw, len(raw)), create_string_buffer(raw, len(raw))
    cn, fn = c_int(-7), c_int(-7)
    clib.strip(cbuf, c_int(len(raw)), byref(cn))
    flib.__m_userio_MOD_strip(fbuf, byref(fn), c_size_t(len(raw)))
    return (cbuf.raw, cn.value), (fbuf.raw, fn.value)


@pytest.mark.parametrize('text', [
    '', ' ', '      ', 'a', ' a', 'a ', '  NACA 0012   ', 'NACA 0012', '   x y  z',
    '\tTAB', 'trailing\t ', ' ' * 47 + 'Z', 'Z' + ' ' * 47,
])
def test_strip_edge_cases(clib, flib, text):
    c, f = run_strip(clib, flib, text)
    assert c == f


def test_strip_random_matches_fortran(clib, flib):
    rng = np.random.default_rng(3)
    alphabet = np.array(list('  ab XYZ09.-'))
    for _ in range(3000):
        n = int(rng.integers(0, 70))
        text = ''.join(rng.choice(alphabet, size=n))
        c, f = run_strip(clib, flib, text)
        assert c == f, repr(text)


def test_strip_result(clib, flib):
    (buf, n), _ = run_strip(clib, flib, '   NACA 0012  ')
    assert buf == b'NACA 0012     ' and n == 9
