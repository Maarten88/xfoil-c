"""Regenerate tests/golden/api_golden.npz from the current xfoil_c build, after checking that every
result the Fortran-backed `xfoil` package can produce is identical (NaN positions equal).

    python tests/golden/make_golden.py        (from xfoil-c/; needs xfoil-python built)

Only regenerate when a clean-up step changes results on purpose (numerics steps), and record why.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
XFOIL_C = HERE.parents[1]
XFOIL_PY = Path(os.environ.get('XFOIL_PYTHON_DIR') or HERE.parents[2] / 'xfoil-python')   # the Fortran-backed original


def same(a, b):
    a, b = np.asarray(a, np.float64), np.asarray(b, np.float64)
    if a.shape != b.shape:
        return False
    d = a.view(np.uint64) != b.view(np.uint64)
    return bool(np.all(np.isnan(a[d]) & np.isnan(b[d])))


def produce(package, out):
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(XFOIL_C), str(XFOIL_PY)]))
    subprocess.run([sys.executable, str(HERE / 'regression_cases.py'), out, '--package', package], check=True, env=env)
    return dict(np.load(out))


if __name__ == '__main__':
    with tempfile.TemporaryDirectory() as td:
        c = produce('xfoil_c', os.path.join(td, 'c.npz'))
        f = produce('xfoil', os.path.join(td, 'f.npz'))
    common = sorted(set(c) & set(f))
    assert set(f) <= set(c), sorted(set(f) - set(c))
    bad = [k for k in common if not same(c[k], f[k])]
    # intentional differences (clean-up C6c): cseq(0.0, 1.2, 0.2) has 6 points with step 0.2 (Fortran package:
    # 5 points with step 0.24, as 1.2/0.2 rounds down), and the cp of that sequence's last point
    INTENTIONAL = {'n23012:003:cseq', 'n23012:004:cp'}
    assert set(bad) <= INTENTIONAL, sorted(set(bad) - INTENTIONAL)
    np.savez_compressed(HERE / 'api_golden.npz', **c)
    nan = sum(int(np.isnan(v).sum()) for v in c.values())
    print('golden: %d arrays (%d also checked against Fortran), %d values, %d NaN'
          % (len(c), len(common), sum(v.size for v in c.values()), nan))
