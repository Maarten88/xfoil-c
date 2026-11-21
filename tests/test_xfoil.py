"""m_xfoil (init) vs the Fortran original: whole-state comparison."""
from ctypes import c_char

import numpy as np
import pytest

from conftest import INDEX_KEYS, _global_buffers, _global_names, globals_diff, to_c_index


def fill_random_state(clib, flib, seed):
    """Identical random contents in every global of both libraries (index globals mapped)."""
    rng = np.random.default_rng(seed)
    clib.xfoil_global_type.restype = c_char
    for k, name in enumerate(_global_names(clib)):
        cbuf, fbuf = _global_buffers(clib, flib, name)
        size = len(fbuf.raw)
        t = clib.xfoil_global_type(k)
        if t in (b'r', b'x'):                    # x: complex = 2 doubles
            f = rng.normal(size=size // 8).astype(np.float64).tobytes()
        elif t == b'i':
            f = rng.integers(1, 60, size=size // 4).astype(np.int32).tobytes()
        elif t == b'l':
            f = rng.integers(0, 2, size=size // 4).astype(np.int32).tobytes()
        else:
            f = bytes(rng.integers(65, 91, size=size).astype(np.uint8))
        if name.endswith(':show_output'):
            f = bytes(size)                 # keep printing off (fixture contract)
        fbuf[:] = f
        cbuf[:] = to_c_index(np.frombuffer(f, np.int32)).tobytes() if name in INDEX_KEYS else f
    assert globals_diff(clib, flib) == []


@pytest.mark.parametrize('seed', range(3))
def test_init_matches_fortran(clib, flib, seed):
    fill_random_state(clib, flib, seed)
    clib.init(); flib.init()
    assert globals_diff(clib, flib) == []


# ------------------------------------------------------------------ naca, pangen

from ctypes import byref, c_int  # noqa: E402

from conftest import global_array, sync_globals  # noqa: E402


def set_both(clib, flib, name, value, dtype=np.float64):
    count = len(_global_buffers(clib, flib, name)[0].raw) // np.dtype(dtype).itemsize
    for lib, fortran in ((clib, False), (flib, True)):
        a = global_array(lib, name, dtype, (count,), fortran)
        a[:len(np.atleast_1d(value))] = value


def init_both(clib, flib, rng):
    sync_globals(clib, flib)
    clib.init(); flib.init()
    set_both(clib, flib, 'npan', int(rng.choice([40, 80, 160, 250, 300])), np.int32)
    set_both(clib, flib, 'cvpar', float(rng.uniform(0.5, 2.0)))
    set_both(clib, flib, 'cterat', float(rng.uniform(0.05, 0.5)))
    set_both(clib, flib, 'ctrrat', float(rng.uniform(0.1, 0.5)))
    if rng.random() < 0.5:                                    # refinement zones
        a, b = sorted(rng.uniform(0.05, 0.95, 2))
        set_both(clib, flib, 'xsref1', a); set_both(clib, flib, 'xsref2', b)
        a, b = sorted(rng.uniform(0.05, 0.95, 2))
        set_both(clib, flib, 'xpref1', a); set_both(clib, flib, 'xpref2', b)
    assert globals_diff(clib, flib) == []


NACA = [12, 2412, 4415, 6, 9, 1408, 23012, 21015, 25012, 24112, 22112, 26012, 30012, 99999]


@pytest.mark.parametrize('k,ides', list(enumerate(NACA)))
def test_naca_matches_fortran(clib, flib, k, ides):
    init_both(clib, flib, np.random.default_rng(k))
    clib.naca(c_int(ides)); flib.__m_xfoil_MOD_naca(byref(c_int(ides)))
    assert globals_diff(clib, flib) == []


def test_naca_nonpositive_returns_without_change(clib, flib):
    """User decision: no prompt; the state stays unchanged (Fortran would ask, so no comparison)."""
    init_both(clib, flib, np.random.default_rng(0))
    clib.naca(c_int(0)); clib.naca(c_int(-12))
    assert globals_diff(clib, flib) == []


def buffer_shape(rng, sharp_le, corners):
    """Biconvex-like buffer airfoil, TE upper -> LE -> TE lower, with optional doubled points."""
    m = int(rng.integers(30, 120))
    t = float(rng.uniform(0.04, 0.15))
    xu = 0.5 * (1 - np.cos(np.linspace(0, np.pi, m)))[::-1]    # 1 -> 0
    if sharp_le:
        yu = 2 * t * xu * (1 - xu)                             # pointed LE (finite slope)
    else:
        yu = t * 5 * (0.2969 * np.sqrt(xu) - 0.126 * xu - 0.3516 * xu**2 + 0.2843 * xu**3 - 0.1015 * xu**4)
    xl, yl = xu[::-1], -yu[::-1] * float(rng.uniform(0.6, 1.0))
    if sharp_le:                                                # doubled LE point
        x, y = np.r_[xu, xl], np.r_[yu, yl]
    else:
        x, y = np.r_[xu, xl[1:]], np.r_[yu, yl[1:]]
    for _ in range(corners):                                    # doubled points elsewhere
        j = int(rng.integers(5, len(x) - 5))
        x, y = np.insert(x, j, x[j]), np.insert(y, j, y[j])
    return x, y


@pytest.mark.parametrize('seed', range(12))
def test_pangen_matches_fortran(clib, flib, seed):
    rng = np.random.default_rng(100 + seed)
    init_both(clib, flib, rng)
    sharp_le, corners = bool(seed % 2), seed % 3
    x, y = buffer_shape(rng, sharp_le, corners)
    set_both(clib, flib, 'xb', x); set_both(clib, flib, 'yb', y)
    set_both(clib, flib, 'nb', len(x), np.int32)
    set_both(clib, flib, 'lbflap', int(seed % 4 == 3), np.int32)
    set_both(clib, flib, 'xbf', 0.7); set_both(clib, flib, 'ybf', 0.0)
    assert globals_diff(clib, flib) == []
    clib.pangen(c_int(1)); flib.__m_xfoil_MOD_pangen(byref(c_int(1)))
    assert globals_diff(clib, flib) == []
    if sharp_le:          # the doubled-LE branch (ible /= 0) was taken
        sb = global_array(clib, 'sb', np.float64, (4 * 640,))[:len(x)]
        sble = global_array(clib, 'sble', np.float64)[0]
        assert np.sum(sb == sble) == 2


def test_pangen_without_buffer(clib, flib):
    init_both(clib, flib, np.random.default_rng(1))
    set_both(clib, flib, 'nb', 1, np.int32)
    clib.pangen(c_int(0)); flib.__m_xfoil_MOD_pangen(byref(c_int(0)))
    assert globals_diff(clib, flib) == []
