# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""Optimization phase: any build variant (XFOIL_ARCH / XFOIL_FAST, see setup.py) against the exact
reference results in tests/golden/api_golden.npz (user criteria, 2026-11-20): small deviations are
fine and individual convergence outcomes may differ, but unconverged points must not become more
frequent. In the exact build this holds trivially (test_regression checks bit-identity there)."""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

GOLDEN = Path(__file__).resolve().parent / 'golden'
REL_TOL = 1e-3          # on values converged in both builds: |a - b| <= REL_TOL * max(|a|, 0.01)


def compare(now, gold):
    """(max scaled deviation over values finite in both, NaN count now, NaN count reference, keys with
    a different shape)."""
    worst, shapes = 0.0, []
    for k, g in gold.items():
        c = now[k]
        if c.shape != g.shape:
            shapes.append(k)
            continue
        m = np.isfinite(c) & np.isfinite(g)
        if m.any():
            worst = max(worst, float(np.max(np.abs(c[m] - g[m]) / np.maximum(np.abs(g[m]), 0.01))))
    nan_now = sum(int(np.isnan(v).sum()) for v in now.values())
    nan_ref = sum(int(np.isnan(v).sum()) for v in gold.values())
    return worst, nan_now, nan_ref, shapes


def test_build_matches_reference_within_tolerance(tmp_path):
    out = tmp_path / 'now.npz'
    env = dict(os.environ, PYTHONPATH=str(GOLDEN.parents[1]))
    subprocess.run([sys.executable, str(GOLDEN / 'regression_cases.py'), str(out)], check=True, env=env)
    now, gold = dict(np.load(out)), dict(np.load(GOLDEN / 'api_golden.npz'))
    assert sorted(now) == sorted(gold)
    worst, nan_now, nan_ref, shapes = compare(now, gold)
    from xfoil_c import build_info
    print('build %s: max scaled deviation %.2e, NaN %d (reference %d)' % (build_info(), worst, nan_now, nan_ref))
    assert not shapes, shapes
    assert worst <= REL_TOL
    assert nan_now <= nan_ref                     # no more unconverged points than the reference
    if build_info() == 'exact':
        assert worst == 0.0 and nan_now == nan_ref
