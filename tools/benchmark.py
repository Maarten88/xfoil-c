"""Wall-clock benchmark of the public API: Fortran-backed `xfoil` vs `xfoil_c` (MIGRATION_PLAN.md C0).

    python -I tools/benchmark.py [--repeat 5]        (from xfoil-c/)

Each package runs in its own fresh process; the best of --repeat runs is reported per workload.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

XFOIL_C = Path(__file__).resolve().parents[1]
# the Fortran-backed original to compare with: $XFOIL_PYTHON_DIR, default xfoil-python next to xfoil-c
XFOIL_PY = Path(os.environ.get('XFOIL_PYTHON_DIR') or XFOIL_C.parent / 'xfoil-python')

WORKER = r'''
import json, sys, time
import numpy as np
sys.path[:0] = [sys.argv[2], sys.argv[3]]
pkg = __import__(sys.argv[1])
from importlib import import_module
Airfoil = import_module(sys.argv[1] + '.model').Airfoil
sys.path.insert(0, sys.argv[4])
from regression_cases import naca0012_coords
repeat = int(sys.argv[5])

def bench(fn):
    best = float('inf')
    for _ in range(repeat):
        t = time.perf_counter(); fn(); best = min(best, time.perf_counter() - t)
    return best

def xf_new():
    xf = pkg.XFoil(); xf.print = False; return xf

def aseq_viscous():
    xf = xf_new(); xf.airfoil = Airfoil(*naca0012_coords()); xf.Re, xf.max_iter = 1e6, 150
    xf.aseq(-20, 20, 0.5)

def cseq_viscous():
    xf = xf_new(); xf.naca('2412'); xf.Re, xf.max_iter = 5e5, 100
    xf.cseq(-0.25, 1.375, 0.0625)    # exact step: same points in both packages (clean-up C6c)

def optimizer_like():
    # many short evaluations on changing airfoils, as an optimizer population does
    xf = xf_new(); xf.Re, xf.max_iter = 7e5, 100
    for spec in ('2410', '2412', '3412', '4412', '2312', '4415', '1408', '23012', '0012', '6409'):
        xf.naca(spec)
        for a in (0.0, 2.0, 4.0, 6.0, 8.0):
            xf.a(a)

def inviscid():
    xf = xf_new(); xf.naca('4415'); xf.Re = 0.0
    for _ in range(20):
        xf.aseq(-5, 15, 1.0)

def repanel():
    xf = xf_new()
    for n in range(60, 300, 10):
        xf.naca('2412'); xf.repanel(n_nodes=n)

print(json.dumps({name: bench(fn) for name, fn in [('aseq viscous (80 pts)', aseq_viscous),
                                                     ('cseq viscous (26 pts)', cseq_viscous),
                                                     ('optimizer-like (10 foils x 5 a)', optimizer_like),
                                                     ('inviscid aseq x20', inviscid),
                                                     ('repanel x24', repanel)]}))
'''


def run(package, repeat):
    p = subprocess.run([sys.executable, '-c', WORKER, package, str(XFOIL_C), str(XFOIL_PY),
                        str(XFOIL_C / 'tests' / 'golden'), str(repeat)],
                       capture_output=True, text=True, check=True, env=dict(os.environ, MPLBACKEND='Agg'))
    return json.loads(p.stdout.strip().splitlines()[-1])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--repeat', type=int, default=5)
    args = ap.parse_args()
    f, c = run('xfoil', args.repeat), run('xfoil_c', args.repeat)
    print('%-34s %12s %12s %8s' % ('workload (best of %d)' % args.repeat, 'Fortran [s]', 'C [s]', 'C/F'))
    for k in f:
        print('%-34s %12.4f %12.4f %8.3f' % (k, f[k], c[k], c[k] / f[k]))
