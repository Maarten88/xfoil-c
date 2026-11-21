# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""Any build variant (XFOIL_ARCH / XFOIL_FAST, see setup.py) and any platform against the exact reference
results in tests/golden/api_golden.npz (user criteria, 2026-11-20): small deviations are fine and individual
convergence outcomes may differ, but unconverged points must not become more frequent.

The results marked unstable in the golden file (regression_cases.py, case 1: a cseq started from a stalled
boundary layer, and the cp after it) are left out here: which of their points converge depends on the last
bits of the arithmetic, so neither tolerance nor NaN count is meaningful for them. In the exact build on the
golden file's platform every result, unstable ones included, must be identical (as test_regression checks)."""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

GOLDEN = Path(__file__).resolve().parent / 'golden'
REL_TOL = 1e-3          # on values converged in both builds: |a - b| <= REL_TOL * max(|a|, 0.01)
sys.path.insert(0, str(GOLDEN))
from regression_cases import platform_tag, results  # noqa: E402


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
    now_all, gold_all = dict(np.load(out)), dict(np.load(GOLDEN / 'api_golden.npz'))
    unstable = set(gold_all['__unstable__'])
    assert set(now_all['__unstable__']) == unstable
    now, gold = results(now_all), results(gold_all)
    assert sorted(now) == sorted(gold)
    from xfoil_c import build_info
    if build_info() == 'exact' and str(gold_all['__platform__']) == platform_tag():
        worst, nan_now, nan_ref, shapes = compare(now, gold)
        assert not shapes and worst == 0.0 and nan_now == nan_ref
        return
    stable = lambda d: {k: v for k, v in d.items() if k not in unstable}      # noqa: E731
    worst, nan_now, nan_ref, shapes = compare(stable(now), stable(gold))
    print('build %s on %s: max scaled deviation %.2e, NaN %d (reference %d; %d unstable results left out)'
          % (build_info(), platform_tag(), worst, nan_now, nan_ref, len(unstable)))
    assert not shapes, shapes
    assert worst <= REL_TOL
    assert nan_now <= nan_ref                     # no more unconverged points than the reference
