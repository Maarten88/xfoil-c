# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""Clean-up C6c (user decisions 2026-11-20): API results that differ from the original package on
purpose."""
import numpy as np
import pytest

import xfoil_c


def fresh(re=0.0):
    xf = xfoil_c.XFoil()
    xf.print = False
    xf.naca('0012')
    xf.Re, xf.max_iter = re, 100
    return xf


@pytest.mark.parametrize('args,expected', [
    ((0, 0.3, 0.1), [0.0, 0.1, 0.2]),                     # original: 0, 0.15
    ((0, 0.7, 0.1), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]), # original: 6 points, step 0.1167
    ((0, 0.75, 0.5), [0.0, 0.5]),                         # original: 0 only (step 0.75)
    ((-2, 4, 1), [-2, -1, 0, 1, 2, 3]),                   # unchanged: range semantics, end excluded
    ((10, 0, -2), [10, 8, 6, 4, 2]),
    ((0, 1, 1), [0.0]),
])
def test_aseq_points_follow_range_semantics(args, expected):
    a = fresh().aseq(*args)[0]
    assert np.allclose(a, expected, atol=1e-12, rtol=0), a


def test_cseq_point_count():
    cl = fresh().cseq(0.0, 1.2, 0.2)[1]                   # 1.2/0.2 rounds to 5.999...: original 5 points
    assert np.allclose(cl, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], atol=1e-9)


def test_cp_distribution_is_nan_without_a_converged_point():
    xf = fresh(re=1e6)
    assert np.isnan(xf.get_cp_distribution()[2]).all()    # no analysis yet
    xf.a(2.0)
    x, y, cp = xf.get_cp_distribution()
    assert np.isfinite(cp).all() and np.isfinite(x).all()
    xf.max_iter = 2
    assert np.isnan(xf.a(18.0)[0])                         # unconverged point
    assert np.isnan(xf.get_cp_distribution()[2]).all()
    xf.max_iter = 100
    xf.a(2.0)
    assert np.isfinite(xf.get_cp_distribution()[2]).all()
    for change in (lambda: setattr(xf, 'M', 0.2), lambda: xf.naca('2412'), lambda: xf.repanel(n_nodes=120)):
        change()                                           # geometry / Mach changed: stale cp
        assert np.isnan(xf.get_cp_distribution()[2]).all()
        xf.a(2.0)
        assert np.isfinite(xf.get_cp_distribution()[2]).all()
    xf.aseq(0, 4, 1)                                       # sequence: validity of its last point
    assert np.isfinite(xf.get_cp_distribution()[2]).all()
