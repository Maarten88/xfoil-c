"""Drop-in check: the Fortran-backed `xfoil` package (xfoil-python source tree) and `xfoil_c` have the
same public surface and give identical results through the public Python API."""
import inspect
import sys
from ctypes import byref, c_bool
from pathlib import Path

import numpy as np
import pytest

import xfoil_c
from conftest import same_reals
from test_api_analysis import naca0012

from conftest import XFOIL_PYTHON, require_xfoil_python  # noqa: E402

require_xfoil_python('xfoil', 'libxfoil_fortran.so', patched=True)          # the reference build
sys.path.insert(0, XFOIL_PYTHON)
import xfoil  # noqa: E402  (loads xfoil-python/xfoil/libxfoil_fortran.so)

from xfoil.model import Airfoil as FAirfoil  # noqa: E402
from xfoil_c.model import Airfoil as CAirfoil  # noqa: E402


EXTENSIONS = ('max_thickness', 'x_max_thickness', 'max_camber', 'x_max_camber', 'le_radius', 'te_angle', 'area',
              'set_hinge', 'flap', 'flap_angle', 'hinge', 'hinge_moment', 'hinge_fx', 'hinge_fy',
              'BL_FIELDS', 'get_bl', 'get_bl_nodes')
# removed from xfoil_c (user decision, 2026-11-20): filter only touched inverse-design state, a no-op
REMOVED = ('filter',)


def public(cls):
    return {n: inspect.signature(getattr(cls, n)) if callable(getattr(cls, n)) else type(getattr(cls, n))
            for n in dir(cls) if not n.startswith('_')}


def test_public_surface_is_identical():
    assert sorted(n for n in dir(xfoil) if not n.startswith('_')) == \
           sorted(n for n in dir(xfoil_c) if not n.startswith('_') and n not in ('lib_version', 'build_info', 'XFoilError'))
    ext = public(xfoil_c.XFoil)
    for name in EXTENSIONS:                         # xfoil_c additions (MIGRATION_PLAN.md Phase 10b)
        ext.pop(name)
    orig = public(xfoil.XFoil)
    for name in REMOVED:
        assert not hasattr(xfoil_c.XFoil, name)
        orig.pop(name)
    assert orig == ext
    assert public(FAirfoil) == public(CAirfoil)
    assert str(inspect.signature(xfoil.XFoil.__init__)) == str(inspect.signature(xfoil_c.XFoil.__init__))


def same(a, b):
    """Equal results: tuples/arrays elementwise, floats bit-identical (NaN == NaN)."""
    if isinstance(a, tuple):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, (FAirfoil, CAirfoil)):
        return same(a.coords, b.coords)
    if isinstance(a, (np.ndarray, float)):
        return same_reals(np.atleast_1d(np.asarray(a, np.float64)), np.atleast_1d(np.asarray(b, np.float64)))
    return a == b


SCRIPTS = {
    'baseline': [('airfoil', 'naca0012'), ('Re', 1e6), ('max_iter', 100), ('a', 10), ('cp_dist',),
                 ('cl', 1.0), ('aseq', -5, 10, 1.0), ('cseq', -0.25, 1.0, 0.25), ('cp_dist',)],   # exact step: count as original (C6c)
    'settings': [('naca', '2412'), ('repanel', dict(n_nodes=200, cv_par=1.2, xt_ref=(0.2, 0.5))), ('Re', 3e5),
                 ('M', 0.3), ('xtr', (0.4, 0.6)), ('n_crit', 7.0), ('max_iter', 60), ('a', 4.0),
                 ('reset_bls',), ('cl', 0.8), ('cp_dist',)],
    'inviscid': [('naca', '4415'), ('Re', 0.0), ('a', 6.0), ('cl', 0.5), ('aseq', 0, 8, 2.0), ('cp_dist',)],
    'unconverged': [('airfoil', 'naca0012'), ('Re', 5e5), ('max_iter', 1), ('aseq', -6, 18, 1.0), ('a', 14.0)],
}


def run(pkg, Airfoil, script):
    xf = pkg.XFoil()
    xf.print = False
    out = []
    for step in script:
        op, args = step[0], step[1:]
        if op == 'airfoil':
            xf.airfoil = Airfoil(*naca0012())
            out.append(xf.airfoil)
        elif op in ('Re', 'M', 'xtr', 'n_crit', 'max_iter'):
            setattr(xf, op, args[0])
            out.append(getattr(xf, op))
        elif op == 'cp_dist':
            out.append(xf.get_cp_distribution())
        elif op == 'repanel':
            out.append(xf.repanel(**args[0]))
            out.append(xf.airfoil)
        else:
            out.append(getattr(xf, op)(*args))
    out.append((xf.Re, xf.M, xf.xtr, xf.n_crit, xf.max_iter, xf.print))
    return out


@pytest.mark.parametrize('name', sorted(SCRIPTS))
def test_public_api_results_identical(clib, flib, name):
    f = run(xfoil, FAirfoil, SCRIPTS[name])
    c = run(xfoil_c, CAirfoil, SCRIPTS[name])
    assert len(f) == len(c)
    for k, (a, b) in enumerate(zip(f, c)):
        assert same(a, b), (name, SCRIPTS[name][k] if k < len(SCRIPTS[name]) else 'getters', a, b)
    for lib in (clib, flib):
        lib.set_print(byref(c_bool(False)))
