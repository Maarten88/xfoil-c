# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""Locate and load the libxfoil_c shared library with ctypes.

Each XFoil instance owns an independent XFOIL state (a C context, `xfoil_t`): `cdll.LoadLibrary`
returns a ContextLibrary that creates a context and passes it as first argument to every API
function (so the calls in xfoil.py are the original ones). The Fortran original shares one state
between all instances.

Hard errors (XFOIL's fatal `stop`, e.g. array limits or duplicated input points) do not end the
process as in the Fortran library: the C call returns early and `cdll` (a ctypes loader whose
functions check the library's error state after each guarded call) raises XFoilError with XFOIL's
message. The library state may be partly updated after such an error; set a new airfoil (or create
a new XFoil) before continuing.
"""
import ctypes
import glob
import os
from ctypes import c_char_p, c_int, create_string_buffer

here = os.path.abspath(os.path.dirname(__file__))

_candidates = [p for pattern in ('libxfoil_c.so', 'libxfoil_c.dylib', 'xfoil_c.dll', 'libxfoil_c.dll')
               for p in glob.glob(os.path.join(here, pattern))]
if not _candidates:
    raise ImportError(
        "libxfoil_c not found in {}. Build it with 'pip install .' or "
        "'python setup.py build_ext --inplace'.".format(here))

lib_path = _candidates[0]
lib = ctypes.cdll.LoadLibrary(lib_path)

lib.xfoil_version.argtypes = []
lib.xfoil_version.restype = c_char_p
lib.xfoil_error.argtypes = [c_char_p, c_int]
lib.xfoil_error.restype = c_int

# public functions that can hit a hard error (run inside the C error guard, src/xerror.h)
GUARDED = frozenset(('init', 'set_airfoil', 'set_naca', 'set_mach', 'repanel', 'alfa', 'cl', 'aseq',
                     'cseq', 'get_cp', 'get_geometry', 'set_hinge', 'flap', 'get_bl'))


class XFoilError(RuntimeError):
    """A hard error in the XFOIL library (a Fortran `stop` in the original), with XFOIL's message."""


def _check(result, func, args):
    buf = create_string_buffer(256)
    if lib.xfoil_error(buf, len(buf)):
        raise XFoilError(buf.value.decode('ascii', 'replace'))
    return result


class CheckedCDLL(ctypes.CDLL):
    """ctypes.CDLL whose GUARDED functions raise XFoilError after a hard library error."""

    def __getitem__(self, name):
        func = super().__getitem__(name)
        if name in GUARDED:
            func.errcheck = _check
        return func


# API functions whose first parameter is the context (src/xfoil_c.h); tests/test_skeleton.py checks
# this list against the header
CONTEXT_FUNCTIONS = frozenset((
    'init', 'set_print', 'get_print', 'set_airfoil', 'get_n_coords', 'get_airfoil', 'set_naca',
    'set_reynolds', 'get_reynolds', 'set_mach', 'get_mach', 'set_xtr', 'get_xtr', 'set_n_crit', 'get_n_crit',
    'set_max_iter', 'get_max_iter', 'reset_bls', 'repanel', 'alfa', 'cl', 'aseq', 'cseq', 'get_n_cp', 'get_cp',
    'get_geometry', 'set_hinge', 'flap', 'get_hinge', 'get_hinge_point', 'clear_hinge', 'get_flap_angle',
    'get_n_bl', 'get_bl'))


class _BoundFunction:
    """A library function called with a context prepended; restype/argtypes/... pass through."""

    def __init__(self, func, ctx):
        object.__setattr__(self, '_func', func)
        object.__setattr__(self, '_ctx', ctx)

    def __call__(self, *args):
        return self._func(self._ctx, *args)

    def __getattr__(self, name):
        return getattr(self._func, name)

    def __setattr__(self, name, value):
        setattr(self._func, name, value)


class ContextLibrary:
    """The library bound to a new, independent XFOIL context (freed with this object)."""

    def __init__(self, path):
        self._dll = CheckedCDLL(path)
        self._name = self._dll._name
        self._dll.xfoil_new.restype = ctypes.c_void_p
        ctx = self._dll.xfoil_new()
        if not ctx:
            raise MemoryError('xfoil_new: cannot allocate an XFOIL context')
        self._ctx = ctypes.c_void_p(ctx)

    def __getattr__(self, name):
        func = getattr(self._dll, name)
        return _BoundFunction(func, self._ctx) if name in CONTEXT_FUNCTIONS else func

    def __del__(self):
        dll, ctx = self.__dict__.get('_dll'), self.__dict__.get('_ctx')
        if dll is not None and ctx is not None:
            dll.xfoil_free(ctx)
            self._ctx = None


class _ContextLoader:
    def LoadLibrary(self, path):
        return ContextLibrary(path)


cdll = _ContextLoader()


def version():
    """str: Version reported by the C library."""
    return lib.xfoil_version().decode('ascii')


lib.xfoil_build_info.argtypes = []
lib.xfoil_build_info.restype = c_char_p


def build_info():
    """str: 'exact' (bit-identical to the Fortran reference) or the build variant (setup.py)."""
    return lib.xfoil_build_info().decode('ascii')
