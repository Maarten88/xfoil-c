# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""Shared fixtures: the C library and, for reference comparisons, the Fortran one.

Fortran module routines are called through their gfortran symbols
(__<module>_MOD_<routine>); every argument is passed by reference.
"""
import os
from ctypes import CDLL, c_bool, c_char, byref

import pytest

# The Fortran reference: an xfoil-python checkout ($XFOIL_PYTHON_DIR, default: next to this repository)
# with the library built (see README, "Tests"). Tests that need it are skipped when it is missing.
XFOIL_PYTHON = os.path.abspath(os.environ.get('XFOIL_PYTHON_DIR') or os.path.join(
    os.path.dirname(__file__), '..', '..', 'xfoil-python'))
# Fortran reference library; override with XFOIL_FORTRAN_LIB=/path/to/libxfoil_fortran.so
FORTRAN_LIB = os.path.abspath(os.environ.get('XFOIL_FORTRAN_LIB') or os.path.join(
    XFOIL_PYTHON, 'xfoil', 'libxfoil_fortran.so'))
TOOLS = os.path.join(os.path.dirname(__file__), '..', 'tools')


def require_xfoil_python(*parts):
    """Path inside the xfoil-python checkout; skips the calling test (module) if it is missing."""
    path = os.path.join(XFOIL_PYTHON, *parts)
    if not os.path.exists(path):
        pytest.skip('xfoil-python reference not found: ' + path, allow_module_level=True)
    return path


def _ctx_functions():
    """Exported routines whose first parameter is the context `xfoil_t *xf` (from the headers)."""
    import re
    src = os.path.join(os.path.dirname(__file__), '..', 'src')
    names = set()
    for f in os.listdir(src):
        if f.endswith('.h'):
            text = open(os.path.join(src, f)).read()
            names |= set(re.findall(r'XFOIL_API\s+[\w\s\*]+?\b(\w+)\s*\(\s*xfoil_t\s*\*\s*xf\b', text))
    return names


class _CtxFunction:
    """A library routine called with the default context prepended (restype/argtypes pass through)."""
    def __init__(self, fn, ctx):
        object.__setattr__(self, '_fn', fn)
        object.__setattr__(self, '_ctx', ctx)

    def __call__(self, *args):
        return self._fn(self._ctx, *args)

    def __getattr__(self, name):
        return getattr(self._fn, name)

    def __setattr__(self, name, value):
        setattr(self._fn, name, value)


class CtxLib:
    """The C library for the tests: internal routines that take `xfoil_t *xf` get the library's
    default context (the one the API functions and xfoil_global use) automatically."""
    def __init__(self, lib):
        from ctypes import c_void_p
        lib.xfoil_default.restype = c_void_p
        self._lib = lib
        self._ctx = c_void_p(lib.xfoil_default())
        self._names = _ctx_functions()

    def __getattr__(self, name):
        fn = getattr(self._lib, name)
        return _CtxFunction(fn, self._ctx) if name in self._names else fn


def build_info():
    from xfoil_c import _lib
    return _lib.build_info()


def pytest_collection_modifyitems(config, items):
    """In a non-exact build (XFOIL_ARCH / XFOIL_FAST, optimization phase) the bit-identical
    comparisons with Fortran and the golden file cannot hold: skip them (every test that uses the
    Fortran library, plus test_regression). tests/test_reference.py checks such builds instead."""
    info = build_info()
    if info == 'exact':
        return
    skip = pytest.mark.skip(reason='bit-identical check; build is "%s" (see test_reference.py)' % info)
    for item in items:
        if 'flib' in getattr(item, 'fixturenames', ()) or item.module.__name__ == 'test_regression':
            item.add_marker(skip)


@pytest.fixture(scope='session')
def clib():
    from xfoil_c import _lib
    lib = CtxLib(_lib.lib)
    lib.get_print.restype = c_bool
    lib.set_print(byref(c_bool(False)))          # on the default context (API functions take it first)
    return lib


@pytest.fixture(scope='session')
def flib():
    if not os.path.exists(FORTRAN_LIB):
        pytest.skip('Fortran reference library not found: ' + FORTRAN_LIB)
    lib = CDLL(FORTRAN_LIB)
    lib.set_print(byref(c_bool(False)))
    # associate i_xbl's pointer aliases (X1 => COM1(1), ...), as Fortran init() does;
    # the C port has them built into the blcom_t union
    lib.__m_xbl_MOD_preptrs()
    return lib


def global_array(lib, name, dtype, shape=None, fortran=False):
    """Writable numpy view of a global (Fortran order), from the C table or a Fortran symbol.

    shape is in Fortran index order, e.g. (IVX, ISX) for ISYS(ivx, isx).
    """
    import numpy as np
    from ctypes import c_char, c_size_t, c_void_p
    if fortran:
        sym = fortran_symbol(name)
        nbytes = np.dtype(dtype).itemsize * (int(np.prod(shape)) if shape else 1)
        buf = (c_char * nbytes).in_dll(lib, sym)
    else:
        lib.xfoil_global.restype = c_void_p
        size = c_size_t()
        addr = lib.xfoil_global(name.encode(), byref(size))
        assert addr, name
        buf = (c_char * size.value).from_address(addr)
    arr = np.frombuffer(buf, dtype=dtype)
    return arr.reshape(shape, order='F') if shape else arr


def _load_index_globals():
    path = os.path.join(TOOLS, 'index_globals.txt')
    with open(path) as f:
        return {l.split('#', 1)[0].strip() for l in f} - {''}


# Globals holding index values (bare name = i_xfoil, else 'module:name'): C stores Fortran k as k-1,
# and Fortran 0 ("none") as -1. INDEX_KEYS has every entry as 'module:name'.
INDEX_GLOBALS = _load_index_globals()
INDEX_KEYS = {k if ':' in k else 'i_xfoil:' + k for k in INDEX_GLOBALS}


def to_c_index(f):
    """Map Fortran 1-based index values (0 = none) to the C convention (0-based, -1 = none)."""
    import numpy as np
    f = np.asarray(f)
    return np.where(f >= 1, f - 1, -1).astype(f.dtype)


def fortran_symbol(name):
    """'module:name' (or bare 'name' = i_xfoil) -> gfortran symbol __<module>_MOD_<name>."""
    module, _, var = name.rpartition(':')
    return '__%s_MOD_%s' % (module or 'i_xfoil', var)


def _global_buffers(clib, flib, name):
    from ctypes import c_char, c_size_t, c_void_p
    clib.xfoil_global.restype = c_void_p
    size = c_size_t()
    addr = clib.xfoil_global(name.encode(), byref(size))
    return (c_char * size.value).from_address(addr), (c_char * size.value).in_dll(flib, fortran_symbol(name))


def _global_names(clib):
    from ctypes import c_char_p
    clib.xfoil_global_name.restype = c_char_p
    return [clib.xfoil_global_name(k).decode() for k in range(clib.xfoil_global_count())]


def sync_globals(clib, flib):
    """Set the whole C i_xfoil state equal to the Fortran state (index globals mapped)."""
    import numpy as np
    for name in _global_names(clib):
        cbuf, fbuf = _global_buffers(clib, flib, name)
        if name in INDEX_KEYS:
            cbuf[:] = to_c_index(np.frombuffer(fbuf.raw, np.int32)).tobytes()
        else:
            cbuf[:] = fbuf.raw


def same_reals(a, b):
    """Doubles (bytes or arrays) bit-identical, except that any NaN equals any NaN in the same place.

    The sign/payload of a NaN is not compared: x86 creates NaNs with the sign bit set and, when two
    NaNs meet, keeps the first operand's, so it depends on the operand order each compiler chooses.
    """
    import numpy as np
    a = np.frombuffer(a, np.float64) if isinstance(a, (bytes, bytearray)) else np.asarray(a, np.float64)
    b = np.frombuffer(b, np.float64) if isinstance(b, (bytes, bytearray)) else np.asarray(b, np.float64)
    if a.shape != b.shape:
        return False
    differ = a.view(np.uint64) != b.view(np.uint64)
    return bool(np.all(np.isnan(a[differ]) & np.isnan(b[differ])))


# Fortran's update places its scratch arrays on VA/VB (EQUIVALENCE) and leaves them there; the C
# update does not (clean-up C4d). So globals_diff skips them; va_vb_diff compares their active
# blocks (system lines 0..nsys-1), which setbl/blsolv define completely.
FORTRAN_SCRATCH = {'i_xfoil:va', 'i_xfoil:vb'}


def va_vb_diff(clib, flib):
    """'va'/'vb' if their blocks 0..nsys-1 differ between C and Fortran (same_reals)."""
    import numpy as np
    nsys = int(global_array(clib, 'nsys', np.int32)[0])
    izx = global_array(clib, 'va', np.float64).size // 6
    out = []
    for name in ('va', 'vb'):
        c = global_array(clib, name, np.float64, (3, 2, izx), False)[:, :, :nsys]
        f = global_array(flib, name, np.float64, (3, 2, izx), True)[:, :, :nsys]
        if not same_reals(np.ascontiguousarray(c), np.ascontiguousarray(f)):
            out.append(name)
    return out


def globals_diff(clib, flib, modules=None):
    """Names ('module:name') of the module globals whose C and Fortran contents differ (index
    globals mapped; real globals compared with same_reals). modules: optional tuple such as
    ('i_xbl', 'i_blpar') to compare only those modules (default: all). VA/VB are skipped
    (FORTRAN_SCRATCH; use va_vb_diff)."""
    import numpy as np
    clib.xfoil_global_type.restype = c_char
    diff = []
    for k, name in enumerate(_global_names(clib)):
        if modules and name.split(':')[0] not in modules:
            continue
        if name in FORTRAN_SCRATCH:
            continue
        cbuf, fbuf = _global_buffers(clib, flib, name)
        cb, fb = cbuf.raw, fbuf.raw
        if name in INDEX_KEYS:
            fb = to_c_index(np.frombuffer(fb, np.int32)).tobytes()
        if cb == fb:
            continue
        if clib.xfoil_global_type(k) in (b'r', b'x') and same_reals(cb, fb):   # x: complex = 2 doubles
            continue
        diff.append(name)
    return diff
