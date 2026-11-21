"""Public-API regression scenarios (clean-up safety net, MIGRATION_PLAN.md C0).

Run in a fresh process (the library state is process-global and `init` does not reset everything):

    python tests/golden/regression_cases.py OUT.npz [--package xfoil_c|xfoil]

Every scenario starts with a new XFoil(), uses only the public Python API and records every returned
value as float64 arrays under "<case>:<n>:<what>". With --package xfoil (the Fortran-backed original,
importable from xfoil-python/) the xfoil_c-only extensions are skipped.
"""
import argparse
import sys
import ast
import importlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                    # naca0012_coords.py


def naca0012_coords():
    """NACA 0012 of the original xfoil-python test suite (embedded copy, naca0012_coords.py)."""
    import naca0012_coords as c
    return [np.array(c.X, dtype=np.float64), np.array(c.Y, dtype=np.float64)]


def random_shape():
    """Deterministic smooth closed shape (TE upper -> LE -> TE lower), PCG64 seed 1."""
    rng = np.random.default_rng(1)
    m = 70
    xu = 0.5 * (1 - np.cos(np.linspace(0, np.pi, m)))[::-1]
    t, c = 0.11, 0.035
    yt = t * 5 * (0.2969 * np.sqrt(xu) - 0.126 * xu - 0.3516 * xu**2 + 0.2843 * xu**3 - 0.1036 * xu**4)
    yc = c * 4 * xu * (1 - xu) * (1 + 0.3 * rng.standard_normal(1)[0] * xu)
    x = np.r_[xu, xu[::-1][1:]]
    y = np.r_[yc + yt, (yc - yt)[::-1][1:]]
    return x, y


class Recorder:
    def __init__(self, xf, ext):
        self.xf, self.ext, self.out, self.k = xf, ext, {}, 0

    def put(self, case, what, value, ext=False):
        """Common results are numbered in call order; xfoil_c-only (ext) results reuse the number of the
        preceding common result with an 'x' suffix, so that the common keys match between packages."""
        if ext:
            key = '%s:%03dx:%s' % (case, self.k - 1, what)
        else:
            key = '%s:%03d:%s' % (case, self.k, what)
            self.k += 1
        assert key not in self.out, key
        self.out[key] = np.atleast_1d(np.asarray(value, dtype=np.float64))

    def state(self, case):
        xf = self.xf
        self.put(case, 'getters', [xf.Re, xf.M, *xf.xtr, xf.n_crit, xf.max_iter, float(xf.print)])
        self.put(case, 'airfoil', xf.airfoil.coords.ravel())
        if self.ext:
            self.put(case, 'geometry', [xf.max_thickness, xf.x_max_thickness, xf.max_camber, xf.x_max_camber,
                                         xf.le_radius, xf.te_angle, xf.area], ext=True)

    def call(self, case, name, *args, **kw):
        res = getattr(self.xf, name)(*args, **kw)
        if res is not None:
            self.put(case, name, np.concatenate([np.atleast_1d(np.asarray(r, dtype=np.float64)) for r in res]))
        if self.ext and name in ('a', 'cl', 'aseq', 'cseq'):
            self.put(case, name + '.hinge', [self.xf.hinge_moment, self.xf.hinge_fx, self.xf.hinge_fy], ext=True)

    def cp(self, case):
        self.put(case, 'cp', np.concatenate(self.xf.get_cp_distribution()))


def run(package):
    pkg = importlib.import_module(package)
    Airfoil = importlib.import_module(package + '.model').Airfoil
    ext = package == 'xfoil_c'
    out = {}

    def new(case):
        xf = pkg.XFoil()
        xf.print = False
        r = Recorder(xf, ext)
        r.case = case
        return xf, r

    # 1: the original test airfoil and calls (xfoil/test.py)
    xf, r = new('naca0012')
    xf.airfoil = Airfoil(*naca0012_coords())
    xf.Re, xf.max_iter = 1e6, 100
    r.state('naca0012'); r.call('naca0012', 'a', 10); r.cp('naca0012'); r.call('naca0012', 'cl', 1.0)
    xf.max_iter = 150
    r.call('naca0012', 'aseq', -20, 20, 0.5); r.call('naca0012', 'cseq', -0.5, 0.5, 0.05); r.cp('naca0012')
    out.update(r.out)

    # 2: cambered NACA, low Re, n_crit and forced transition
    xf, r = new('n2412')
    xf.naca('2412'); xf.Re, xf.n_crit, xf.xtr, xf.max_iter = 3e5, 7.0, (0.3, 1.0), 100
    r.state('n2412'); r.call('n2412', 'aseq', -5, 15, 1.0); r.call('n2412', 'cl', 1.2); r.cp('n2412')
    out.update(r.out)

    # 3: inviscid
    xf, r = new('n4415')
    xf.naca('4415'); xf.Re = 0.0
    r.state('n4415'); r.call('n4415', 'a', 5.0); r.call('n4415', 'cl', 0.8); r.call('n4415', 'aseq', 0, 10, 2.0)
    r.cp('n4415')
    out.update(r.out)

    # 4: compressible, 5-digit
    xf, r = new('n23012')
    xf.naca('23012'); xf.Re, xf.M, xf.max_iter = 2e6, 0.3, 100
    r.state('n23012'); r.call('n23012', 'aseq', -4, 12, 2.0); r.call('n23012', 'cseq', 0.0, 1.2, 0.2)
    r.cp('n23012')
    out.update(r.out)

    # 5: repaneling options, reset_bls
    xf, r = new('repanel')
    xf.naca('2412'); xf.Re, xf.max_iter = 5e5, 100
    xf.repanel(n_nodes=240, cv_par=1.5, cte_ratio=0.25, ctr_ratio=0.3, xt_ref=(0.2, 0.6), xb_ref=(0.3, 0.7))
    r.state('repanel'); r.call('repanel', 'aseq', 0, 8, 1.0)
    xf.reset_bls(); r.call('repanel', 'a', 3.0)
    r.call('repanel', 'a', 3.0); r.cp('repanel')
    xf.repanel(); r.state('repanel'); r.call('repanel', 'a', 3.0)
    out.update(r.out)

    # 6: unconverged points (NaN results)
    xf, r = new('unconv')
    xf.airfoil = Airfoil(*naca0012_coords()); xf.Re, xf.max_iter = 5e5, 1
    r.call('unconv', 'aseq', -6, 18, 1.0); r.call('unconv', 'a', 14.0); r.call('unconv', 'cl', 1.5)
    out.update(r.out)

    # 7: a non-NACA shape, clockwise input, repaneled
    xf, r = new('shape')
    x, y = random_shape()
    xf.airfoil = Airfoil(x[::-1].copy(), y[::-1].copy())
    xf.repanel(n_nodes=180); xf.Re, xf.max_iter = 8e5, 100
    r.state('shape'); r.call('shape', 'aseq', -2, 10, 1.0); r.call('shape', 'cl', 0.9); r.cp('shape')
    out.update(r.out)

    # 8: extensions: hinge moment and flap (xfoil_c only)
    if ext:
        xf, r = new('flap')
        xf.naca('2412'); xf.Re, xf.max_iter = 1e6, 100
        xf.set_hinge(0.75, 0.5, y_relative=True)
        r.put('flap', 'hinge', xf.hinge); r.call('flap', 'a', 4.0)
        xf.flap(5.0); r.state('flap'); r.call('flap', 'a', 4.0); r.cp('flap')
        xf.flap(-12.0); r.call('flap', 'aseq', 0, 6, 2.0)      # absolute: -12 deg (was cumulative -7 before 2026-11-21)
        xf.Re = 0.0; r.call('flap', 'cl', 0.5)
        out.update(r.out)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--package', default='xfoil_c', choices=('xfoil_c', 'xfoil'))
    args = ap.parse_args()
    np.savez(args.out, **run(args.package))
