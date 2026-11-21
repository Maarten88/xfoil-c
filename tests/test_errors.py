"""Clean-up C2: hard library errors (Fortran `stop`) raise xfoil_c.XFoilError with XFOIL's message instead of
ending the process; the library is usable again afterwards."""
import numpy as np
import pytest

import xfoil_c
from xfoil_c.model import Airfoil
from test_api_analysis import naca0012


def fresh():
    xf = xfoil_c.XFoil()
    xf.print = False
    return xf


def baseline_a10():
    xf = fresh()
    xf.airfoil = Airfoil(*naca0012())
    xf.Re, xf.max_iter = 1e6, 100
    return xf.a(10)


def test_duplicated_first_point_raises():
    x, y = naca0012()
    xf = fresh()
    with pytest.raises(xfoil_c.XFoilError, match='SEGSPL:  First input point duplicated'):
        xf.airfoil = Airfoil(np.r_[x[0], x], np.r_[y[0], y])


def test_too_many_input_points_raises():
    t = np.linspace(0, np.pi, 600)
    xu = 0.5 * (1 + np.cos(t))
    yu = 0.06 * np.sin(t) * (1 - xu)
    x, y = np.r_[xu, xu[::-1][1:]], np.r_[yu, -yu[::-1][1:]]        # 1199 points > NMAX = 1000
    xf = fresh()
    with pytest.raises(xfoil_c.XFoilError, match='SPLIND: array overflow'):
        xf.airfoil = Airfoil(x, y)


def test_too_many_panels_for_the_solver_raises():
    xf = fresh()
    xf.naca('0012')
    xf.repanel(n_nodes=600)                                         # > NVX = 500 unknowns in ludcmp
    with pytest.raises(xfoil_c.XFoilError, match='LUDCMP: Array overflow'):
        xf.a(2.0)


def test_library_usable_after_errors():
    expected = baseline_a10()
    test_duplicated_first_point_raises()
    test_too_many_panels_for_the_solver_raises()
    assert baseline_a10() == expected
    assert np.round(expected[:4], 4).tolist() == [1.0809, 0.0150, 0.0053, -5.5077]


def test_error_is_cleared_after_reading():
    xf = fresh()
    with pytest.raises(xfoil_c.XFoilError):
        xf.airfoil = Airfoil(np.array([0.0, 0.0, 1.0, 0.0]), np.array([0.0, 0.0, 0.1, 0.0]))
    xf.naca('2412')                                                 # no stale error raised
    assert xf.airfoil.n_coords > 100
