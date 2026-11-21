# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""m_xsolve.f90 ludcmp/baksub -> src/xsolve.c: bit-identical to Fortran (C pivots are 0-based)."""
import subprocess
import sys
import textwrap
from ctypes import byref, c_int, c_void_p

import numpy as np

from conftest import FORTRAN_LIB


def ptr(arr):
    return arr.ctypes.data_as(c_void_p)


def factor_and_solve(lib, fortran, a, b, nsiz):
    n = a.shape[0]
    A = np.zeros((nsiz, nsiz), order='F')
    A[:n, :n] = a
    indx = np.zeros(nsiz, dtype=np.int32)
    B = np.zeros(nsiz)
    B[:n] = b
    if fortran:
        lib.__m_xsolve_MOD_ludcmp(byref(c_int(nsiz)), byref(c_int(n)), ptr(A), ptr(indx))
        lib.__m_xsolve_MOD_baksub(byref(c_int(nsiz)), byref(c_int(n)), ptr(A), ptr(indx), ptr(B))
        indx[:n] -= 1
    else:
        lib.ludcmp(c_int(nsiz), c_int(n), ptr(A), ptr(indx))
        lib.baksub(c_int(nsiz), c_int(n), ptr(A), ptr(indx), ptr(B))
    return A.tobytes(), indx[:n].tobytes(), B.tobytes()


def test_ludcmp_baksub_match_fortran(clib, flib):
    pivots = leading_zero_rhs = 0
    for seed in range(120):
        rng = np.random.default_rng(seed)
        n = int(rng.integers(1, 120))
        nsiz = n + int(rng.integers(0, 5))
        a = rng.normal(size=(n, n))
        if seed % 3 == 0:
            a += np.diag(np.full(n, 3.0 * n))         # diagonally dominant: few pivots
        b = rng.normal(size=n)
        if seed % 4 == 0:
            b[: n // 2] = 0.0                          # exercises baksub's ii logic
            leading_zero_rhs += 1
        c = factor_and_solve(clib, False, a, b, nsiz)
        assert c == factor_and_solve(flib, True, a, b, nsiz), seed
        pivots += int(np.sum(np.frombuffer(c[1], np.int32) != np.arange(n)))
        x = np.frombuffer(c[2])[:n]
        assert np.allclose(a @ x, b, atol=1e-6 * (1 + np.abs(b).max()))
    assert pivots > 1000 and leading_zero_rhs > 0


def test_ludcmp_large(clib, flib):
    """n = 500 (the NVX limit) in an nsiz = IQX = 640 array, as used for AIJ."""
    n = 500
    rng = np.random.default_rng(n)
    a, b = rng.normal(size=(n, n)), rng.normal(size=n)
    assert factor_and_solve(clib, False, a, b, 640) == factor_and_solve(flib, True, a, b, 640)


STOP_SCRIPT = textwrap.dedent('''
    import sys
    from ctypes import CDLL, byref, c_int, c_void_p
    import numpy as np
    lib = CDLL(sys.argv[1]); a = np.zeros((501, 501), order='F'); ix = np.zeros(501, dtype=np.int32)
    P = lambda v: v.ctypes.data_as(c_void_p)
    if sys.argv[2] == 'c':
        lib.ludcmp(c_int(501), c_int(501), P(a), P(ix))
    else:
        lib.__m_xsolve_MOD_ludcmp(byref(c_int(501)), byref(c_int(501)), P(a), P(ix))
    print('NOT STOPPED')
''')


def test_ludcmp_stop(clib, flib):
    from xfoil_c import _lib
    res = {}
    for tag, path in (('c', _lib.lib_path), ('f', FORTRAN_LIB)):
        p = subprocess.run([sys.executable, '-c', STOP_SCRIPT, path, tag], capture_output=True, text=True)
        res[tag] = (p.returncode, p.stdout.strip(), p.stderr.strip())
    assert res['c'] == res['f'] == (0, '', 'STOP LUDCMP: Array overflow. Increase NVX.')


def test_gauss_matches_fortran(clib, flib):
    for seed in range(300):
        rng = np.random.default_rng(seed)
        nn = int(rng.integers(1, 12))
        nsiz = nn + int(rng.integers(0, 3))
        nrhs = int(rng.integers(1, 4))
        z = np.asfortranarray(rng.normal(size=(nsiz, nsiz)))
        r = np.asfortranarray(rng.normal(size=(nsiz, nrhs)))
        out = []
        for fortran in (False, True):
            zz, rr = z.copy(order='F'), r.copy(order='F')
            if fortran:
                flib.__m_xsolve_MOD_gauss(byref(c_int(nsiz)), byref(c_int(nn)), ptr(zz), ptr(rr), byref(c_int(nrhs)))
            else:
                clib.gauss(c_int(nsiz), c_int(nn), ptr(zz), ptr(rr), c_int(nrhs))
            out.append(zz.tobytes() + rr.tobytes())
        assert out[0] == out[1], seed


def test_blsolv_matches_fortran(clib, flib):
    """Random diagonally dominant block systems; isys/iblte written 1-based (F) / 0-based (C)."""
    from conftest import global_array, globals_diff, sync_globals, to_c_index, va_vb_diff
    IZX, IVX = 722, 452
    sync_globals(clib, flib)
    branch = {'small': 0, 'large': 0}
    for seed in range(40):
        rng = np.random.default_rng(100 + seed)
        nsys = int(rng.integers(6, 120))
        n_te1 = int(rng.integers(2, nsys - 3))                  # side-1 lines: 1..n_te1
        isys = np.zeros((IVX, 2), dtype=np.int32, order='F')
        isys[1:n_te1 + 1, 0] = np.arange(1, n_te1 + 1)           # stations 2..n_te1+1 on side 1
        isys[1:nsys - n_te1 + 1, 1] = np.arange(n_te1 + 1, nsys + 1)
        iblte = np.array([n_te1 + 1, int(rng.integers(2, nsys - n_te1))], dtype=np.int32)
        va = np.zeros((3, 2, IZX), order='F'); vb = np.zeros((3, 2, IZX), order='F')
        vdel = np.zeros((3, 2, IZX), order='F'); vm = np.zeros((3, IZX, IZX), order='F')
        va[:, :, :nsys] = rng.normal(size=(3, 2, nsys)) * 0.3
        va[0, 0, :nsys] += 4.0
        va[1, 1, :nsys] += 4.0
        vb[:, :, :nsys] = rng.normal(size=(3, 2, nsys)) * 0.3
        vm[:, :nsys, :nsys] = rng.normal(size=(3, nsys, nsys)) * rng.choice([1e-6, 1e-3, 0.05], size=(3, nsys, nsys))
        for iv in range(nsys):
            vm[2, iv, iv] += 4.0
        vdel[:, :, :nsys] = rng.normal(size=(3, 2, nsys))
        vz = rng.normal(size=(3, 2)) * 0.3
        vaccel = float(rng.choice([1e-4, 1e-2]))
        for lib, fortran in ((clib, False), (flib, True)):
            for name, arr in (('va', va), ('vb', vb), ('vm', vm), ('vdel', vdel), ('vz', vz)):
                global_array(lib, name, np.float64, arr.shape, fortran)[...] = arr
            global_array(lib, 'nsys', np.int32, None, fortran)[:] = nsys
            global_array(lib, 'n', np.int32, None, fortran)[:] = 10
            global_array(lib, 's', np.float64, (IZX,), fortran)[:10] = np.linspace(0.0, 2.0, 10)
            global_array(lib, 'vaccel', np.float64, None, fortran)[:] = vaccel
            global_array(lib, 'isys', np.int32, (IVX, 2), fortran)[...] = isys if fortran else to_c_index(isys)
            global_array(lib, 'iblte', np.int32, (2,), fortran)[:] = iblte if fortran else iblte - 1
        assert globals_diff(clib, flib) == []
        clib.blsolv()
        flib.__m_xsolve_MOD_blsolv()
        assert globals_diff(clib, flib) == [], seed
        assert va_vb_diff(clib, flib) == [], seed
        small = np.abs(vm[:, :nsys, :nsys]) <= vaccel
        branch['small'] += int(small.sum() > 0)
        branch['large'] += int((~small).sum() > 0)
    assert min(branch.values()) > 10
