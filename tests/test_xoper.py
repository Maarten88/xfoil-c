"""m_xoper (viscal) and s_xoper (mhinge) vs the Fortran originals: whole-state comparison."""
from ctypes import byref, c_int

import numpy as np
import pytest

from conftest import global_array, globals_diff
import test_xbl as XB


def prepare_viscal(clib, flib, seed, closed_te, sharp_le):
    """prepare_setbl state + the flags and constants viscal reads; branches varied by seed."""
    st, n = XB.prepare_setbl(clib, flib, seed, closed_te, sharp_le)
    rng = np.random.default_rng(60000 + seed)
    st.set('dtor', np.pi / 180.0)
    st.set('minf_cl', 0.0)
    st.set('reinf_cl', float(global_array(clib, 'reinf1', np.float64)[0]))
    st.set('reinf', float(global_array(clib, 'reinf1', np.float64)[0]))
    st.set('minf', float(global_array(clib, 'minf1', np.float64)[0]))
    st.set('xcmref', 0.25)
    st.set('ycmref', 0.0)
    flags = {'lwake': 1, 'lwdij': 1, 'ladij': 1, 'lvisc': 1,
             'lipan': int(seed % 3 != 0),          # 0: redo stfind/iblpan/xicalc/iblsys
             'lvconv': int(seed % 4 == 1),         # 1: CL/CD from the current state first
             'lflap': int(seed % 5 == 2)}          # 1: hinge moment (mhinge)
    for name, v in flags.items():
        st.set(name, v, np.int32)
    x, _ = st.views('x')
    st.set('xof', float(np.min(x[:n]) + rng.uniform(0.6, 0.9) * np.ptp(x[:n])))
    st.set('yof', float(rng.uniform(-0.01, 0.01)))
    assert globals_diff(clib, flib) == []
    return st, n


@pytest.mark.parametrize('seed,closed_te,sharp_le', XB.MARCH_STATES)
def test_viscal_matches_fortran(clib, flib, seed, closed_te, sharp_le):
    prepare_viscal(clib, flib, seed, closed_te, sharp_le)
    niter = [1, 5, 20][seed % 3]
    flib.__m_xoper_MOD_viscal.restype = c_int
    ok_f = flib.__m_xoper_MOD_viscal(byref(c_int(niter)))
    ok_c = clib.viscal(c_int(niter))
    assert ok_c == ok_f
    assert XB.diff_allowing_nan_states(clib, flib) == []


def test_viscal_reaches_branches(clib, flib):
    """Some states converge (lvconv set by viscal) and some abort or fail."""
    conv = []
    for seed, cte, sle in XB.MARCH_STATES:
        prepare_viscal(clib, flib, seed, cte, sle)
        global_array(clib, 'lvconv', np.int32)[:] = 0
        ok = clib.viscal(c_int(30))
        conv.append((ok, int(global_array(clib, 'lvconv', np.int32)[0])))
    assert (1, 1) in conv and any(c != (1, 1) for c in conv), conv


@pytest.mark.parametrize('seed,closed_te,sharp_le', XB.MARCH_STATES[:8])
@pytest.mark.parametrize('lvisc', [0, 1])
def test_mhinge_matches_fortran(clib, flib, seed, closed_te, sharp_le, lvisc):
    st, n = prepare_viscal(clib, flib, seed, closed_te, sharp_le)
    rng = np.random.default_rng(70000 + seed)
    st.set('lflap', 1, np.int32)
    st.set('lvisc', lvisc, np.int32)
    st.set('cpi', rng.normal(size=st.views('cpi')[0].size))
    st.set('cpv', rng.normal(size=st.views('cpv')[0].size))
    clib.mhinge(); flib.__s_xoper_MOD_mhinge()
    assert globals_diff(clib, flib) == []
    assert global_array(clib, 'hmom', np.float64)[0] != 0.0


# ------------------------------------------------------------------ specal, speccl, fcpmin

def prepare_inviscid(clib, flib, seed, closed_te, sharp_le):
    st, n = prepare_viscal(clib, flib, seed, closed_te, sharp_le)
    rng = np.random.default_rng(80000 + seed)
    for name, v in (('lgamu', int(seed % 2)), ('lqaij', int(seed % 2)),   # 0: ggcalc recomputes
                    ('matyp', 1 + seed % 3), ('retyp', 1 + (seed // 3) % 3),
                    ('lvisc', int(seed % 4 < 2)), ('lflap', int(seed % 5 == 2))):
        st.set(name, v, np.int32)
    st.set('minf1', float(rng.choice([0.0, rng.uniform(0.1, 0.75)])))
    st.set('clspec', float(rng.uniform(-0.3, 1.5)))
    st.set('alfa', float(rng.uniform(-0.15, 0.25)))
    assert globals_diff(clib, flib) == []
    return st, n


@pytest.mark.parametrize('seed,closed_te,sharp_le', XB.MARCH_STATES)
@pytest.mark.parametrize('routine', ['specal', 'speccl'])
def test_specal_speccl_match_fortran(clib, flib, seed, closed_te, sharp_le, routine):
    prepare_inviscid(clib, flib, seed, closed_te, sharp_le)
    getattr(clib, routine)(); getattr(flib, '__m_xoper_MOD_' + routine)()
    assert globals_diff(clib, flib) == []


def test_specal_reaches_mach_relaxation(clib, flib):
    """Some states halve rlx in specal's M(CL) loop (matyp > 1, Mach limit hit)."""
    rlx = []
    for seed, cte, sle in XB.MARCH_STATES:
        prepare_inviscid(clib, flib, seed, cte, sle)
        clib.specal()
        rlx.append(float(global_array(clib, 'rlx', np.float64)[0]))
    assert min(rlx) < 1.0 and max(rlx) == 1.0, rlx


@pytest.mark.parametrize('seed', range(6))
def test_fcpmin_matches_fortran(clib, flib, seed):
    st, _ = prepare_viscal(clib, flib, seed, False, False)
    rng = np.random.default_rng(90000 + seed)
    st.set('lvisc', seed % 2, np.int32)
    st.set('cpi', rng.normal(size=st.views('cpi')[0].size))
    st.set('cpv', rng.normal(size=st.views('cpv')[0].size))
    clib.fcpmin(); flib.__m_xoper_MOD_fcpmin()
    assert globals_diff(clib, flib) == []
