# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""m_xblsys.f90 -> src/xblsys.c. Part 1: closure correlations, swept across their branch
boundaries; every output must be bit-identical to Fortran (NaN rule: conftest.same_reals)."""
import itertools
from ctypes import byref, c_double, c_int

import numpy as np
import pytest

from conftest import global_array, same_reals

# routine -> number of scalar inputs (by value in C, by reference in Fortran) and outputs
SIGS = {'hkin': (2, 3), 'dil': (2, 3), 'dilw': (2, 3), 'hsl': (3, 4), 'cfl': (3, 4),
        'hst': (3, 4), 'cft': (3, 4), 'hct': (2, 3)}

HK = np.concatenate([np.linspace(1.02, 12.0, 160), [4.0, 4.35, 5.5, 3.5, 4.0 - 1e-12, 4.35 + 1e-12, 5.5 - 1e-12,
                                                     16.0, 20.0]])
RT = np.array([1.0, 10.0, 150.0, 200.0, 200.0 + 1e-9, 300.0, 400.0, 400.0 + 1e-9, 1e3, 1e4, 1e6])
MSQ = np.array([0.0, 0.05, 0.3, 0.8, 1.2])


def call(lib, fortran, name, args):
    nin, nout = SIGS[name]
    outs = [c_double(-1.0) for _ in range(nout)]
    if fortran:
        getattr(lib, '__m_xblsys_MOD_' + name)(*[byref(c_double(a)) for a in args], *[byref(o) for o in outs])
    else:
        getattr(lib, name)(*[c_double(a) for a in args], *[byref(o) for o in outs])
    return [o.value for o in outs]


def inputs(name):
    if name == 'hkin':                                   # (h, msq)
        return itertools.product(np.linspace(1.0, 8.0, 80), MSQ)
    if name == 'hct':                                    # (hk, msq)
        return itertools.product(HK, MSQ)
    if name in ('dil', 'dilw'):                          # (hk, rt)
        return itertools.product(HK, RT)
    return itertools.product(HK, RT, MSQ)               # (hk, rt, msq)


@pytest.fixture
def cffac_set(clib, flib):
    for lib, fortran in ((clib, False), (flib, True)):
        global_array(lib, 'i_blpar:cffac', np.float64, None, fortran)[:] = 1.0
    return 1.0


@pytest.mark.parametrize('name', sorted(SIGS))
def test_closure_matches_fortran(clib, flib, cffac_set, name):
    n = 0
    for args in inputs(name):
        c = call(clib, False, name, args)
        f = call(flib, True, name, args)
        assert same_reals(np.array(c), np.array(f)), (name, args, c, f)
        n += 1
    assert n > 300


def test_branches_reached(clib, flib, cffac_set):
    """hst: attached/separated and both Rt limiters; cft: grt and arg limiters."""
    hk_ho = {'attached': 0, 'separated': 0}
    for hk, rt, msq in inputs('hst'):
        ho = 3.0 + 400.0 / rt if rt > 400.0 else 4.0
        hk_ho['attached' if hk < ho else 'separated'] += 1
    assert min(hk_ho.values()) > 100
    assert any(np.log(rt) < 3.0 for rt in RT)                      # cft: grt clamped to 3
    assert any(-1.33 * hk < -20.0 for hk in HK)                     # cft: arg clamped to -20


# ------------------------------------------------------------------ part 2: amplification

def dampl_call(lib, fortran, name, hk, th, rt):
    outs = [c_double(-1.0) for _ in range(4)]
    if fortran:
        getattr(lib, '__m_xblsys_MOD_' + name)(byref(c_double(hk)), byref(c_double(th)), byref(c_double(rt)),
                                               *[byref(o) for o in outs])
    else:
        getattr(lib, name)(c_double(hk), c_double(th), c_double(rt), *[byref(o) for o in outs])
    return np.array([o.value for o in outs])


DAMPL_HK = np.concatenate([np.linspace(1.05, 20.0, 120), np.linspace(3.3, 4.2, 40)])
DAMPL_TH = np.array([1e-4, 1e-3, 1e-2])


def grcrit(hk):
    hmi = 1.0 / (hk - 1.0)
    return 2.492 * hmi ** 0.43 + 0.7 * (np.tanh(14.0 * hmi - 9.24) + 1.0)


@pytest.mark.parametrize('name', ['dampl', 'dampl2'])
def test_dampl_matches_fortran(clib, flib, name):
    counts = {'off': 0, 'ramp': 0, 'full': 0, 'blend': 0}
    for hk in DAMPL_HK:
        g = grcrit(hk)
        # Rtheta: far below, inside the +/- dgr ramp, and above Rcrit
        for gr in np.concatenate([[g - 0.5, g - 0.079, g + 0.079, g + 0.5, g + 2.0], np.linspace(g - 0.08, g + 0.08, 7)]):
            rt = 10.0 ** gr
            for th in DAMPL_TH:
                c = dampl_call(clib, False, name, hk, th, rt)
                f = dampl_call(flib, True, name, hk, th, rt)
                assert same_reals(c, f), (name, hk, th, rt, c, f)
            rn = (np.log10(rt) - (g - 0.08)) / 0.16
            counts['off' if rn < 0 else 'ramp' if rn < 1 else 'full'] += 1
            counts['blend'] += 3.5 <= hk < 4.0
    assert min(counts.values()) > 20, counts


def axset_call(lib, fortran, args, idampv):
    outs = [c_double(-1.0) for _ in range(9)]
    if fortran:
        lib.__m_xblsys_MOD_axset(*[byref(c_double(a)) for a in args], byref(c_int(idampv)), *[byref(o) for o in outs])
    else:
        lib.axset(*[c_double(a) for a in args], c_int(idampv), *[byref(o) for o in outs])
    return np.array([o.value for o in outs])


@pytest.mark.parametrize('idampv', [0, 1])
def test_axset_matches_fortran(clib, flib, idampv):
    rng = np.random.default_rng(idampv)
    near = 0
    for _ in range(3000):
        hk1, hk2 = rng.uniform(1.05, 8.0, 2)
        t1, t2 = 10.0 ** rng.uniform(-4, -2, 2)
        rt1, rt2 = 10.0 ** rng.uniform(1, 5, 2)
        acrit = 9.0
        a1, a2 = rng.uniform(0.0, 12.0, 2) if rng.random() < 0.5 else rng.uniform(8.0, 10.0, 2)
        args = (hk1, t1, rt1, a1, hk2, t2, rt2, a2, acrit)
        assert same_reals(axset_call(clib, False, args, idampv), axset_call(flib, True, args, idampv)), args
        near += 0 < 20.0 * (acrit - 0.5 * (a1 + a2)) < 20.0
    assert near > 300


# ------------------------------------------------------------------ part 3: station variables

from conftest import globals_diff, sync_globals  # noqa: E402

BL = ('i_xbl', 'i_blpar')


def set_bl_parameters(clib, flib, rng):
    """BL parameter defaults (blpini) and compressibility parameters for a random Mach/Re."""
    gacon, gbcon = 6.70, 0.75
    minf = float(rng.choice([0.0, rng.uniform(0.05, 0.7)]))
    reinf = 10.0 ** rng.uniform(4.5, 7)
    gm1 = 0.4
    beta = np.sqrt(1.0 - minf ** 2)
    hstinv = gm1 * minf ** 2 / (1.0 + 0.5 * gm1 * minf ** 2)
    herat = 1.0 - 0.5 * hstinv
    values = {
        'i_blpar:sccon': 5.6, 'i_blpar:gacon': gacon, 'i_blpar:gbcon': gbcon, 'i_blpar:gccon': 18.0,
        'i_blpar:dlcon': 0.9, 'i_blpar:ctrcon': 1.8, 'i_blpar:ctrcex': 3.3, 'i_blpar:duxcon': 1.0,
        'i_blpar:ctcon': 0.5 / (gacon ** 2 * gbcon), 'i_blpar:cffac': 1.0,
        'i_xbl:gm1bl': gm1, 'i_xbl:qinfbl': 1.0,
        'i_xbl:tkbl': minf ** 2 / (1.0 + beta) ** 2, 'i_xbl:tkbl_ms': float(rng.uniform(0.1, 0.5)),
        'i_xbl:hstinv': hstinv, 'i_xbl:hstinv_ms': float(rng.uniform(0.1, 0.5)),
        'i_xbl:rstbl': (1.0 + 0.5 * gm1 * minf ** 2) ** (1.0 / gm1), 'i_xbl:rstbl_ms': float(rng.uniform(0.0, 0.5)),
        'i_xbl:hvrat': 0.35,
        'i_xbl:reybl': reinf * np.sqrt(herat ** 3) * 1.35 / (herat + 0.35),
        'i_xbl:reybl_re': float(rng.uniform(0.5, 1.0)), 'i_xbl:reybl_ms': float(rng.uniform(-1.0, 1.0)),
    }
    for name, val in values.items():
        for lib, fortran in ((clib, False), (flib, True)):
            global_array(lib, name, np.float64, None, fortran)[:] = val


def blprv_call(clib, flib, args):
    clib.blprv(*[c_double(a) for a in args])
    flib.__m_xblsys_MOD_blprv(*[byref(c_double(a)) for a in args])


def station(rng, ityp):
    """Primary variables (xsi, ampl, ctau, theta, dstar, dswake, uei) over a wide range."""
    th = 10.0 ** rng.uniform(-5, -2)
    # H near 1 reaches the hk >= 1.05 / 1.00005 floors and the Us clamp (H* -> 2 as hk -> 1)
    hk = float(rng.choice([rng.uniform(1.0, 1.15), rng.uniform(1.15, 3.0), rng.uniform(3.0, 8.0),
                           rng.uniform(8.0, 30.0)]))
    dw = float(rng.uniform(0, 0.5) * th) if ityp == 3 else 0.0
    return (float(rng.uniform(0.0, 1.0)), float(rng.uniform(0.0, 12.0)), float(rng.uniform(0.0, 0.4)),
            th, hk * th + dw, dw, float(rng.uniform(0.05, 1.6)))


def set_logical(clib, flib, name, value):
    for lib, fortran in ((clib, False), (flib, True)):
        global_array(lib, name, np.int32, None, fortran)[:] = value


def blcom_fields():
    """blcom_t field -> 0-based slot, read from the generated src/ixbl.h."""
    import re
    from pathlib import Path
    hdr = (Path(__file__).resolve().parents[1] / 'src' / 'ixbl.h').read_text()
    return {f: int(k) - 1 for f, k in re.findall(r'double (\w+);\s*/\* \w+ / \w+ = COM\((\d+)\)', hdr)}


F = blcom_fields()


def test_station_variables_match_fortran(clib, flib):
    sync_globals(clib, flib)
    hits = {'us_clamp': 0, 'de_clamp': 0, 'hk_min': 0}
    for case in range(150):
        rng = np.random.default_rng(5000 + case)
        set_bl_parameters(clib, flib, rng)
        for ityp in (1, 2, 3):
            # station "1": primary -> secondary, then COM1 = COM2 (as the marching routines do)
            blprv_call(clib, flib, station(rng, ityp))
            clib.blkin(); flib.__m_xblsys_MOD_blkin()
            clib.blvar(c_int(ityp)); flib.__m_xblsys_MOD_blvar(byref(c_int(ityp)))
            for lib, fortran in ((clib, False), (flib, True)):
                c1 = global_array(lib, 'i_xbl:com1', np.float64, (73,), fortran)
                c2 = global_array(lib, 'i_xbl:com2', np.float64, (73,), fortran)
                c1[:] = c2
            assert globals_diff(clib, flib, BL) == [], ('station 1', case, ityp)
            # station "2", each step compared
            args = station(rng, ityp)
            blprv_call(clib, flib, args)
            assert globals_diff(clib, flib, BL) == [], ('blprv', case, ityp)
            clib.blkin(); flib.__m_xblsys_MOD_blkin()
            assert globals_diff(clib, flib, BL) == [], ('blkin', case, ityp)
            clib.blvar(c_int(ityp)); flib.__m_xblsys_MOD_blvar(byref(c_int(ityp)))
            assert globals_diff(clib, flib, BL) == [], ('blvar', case, ityp)
            c2 = global_array(clib, 'i_xbl:com2', np.float64, (73,), False)
            hits['us_clamp'] += bool(c2[F['us']] in (0.98, 0.99995))
            hits['de_clamp'] += bool(c2[F['de']] == 12.0 * c2[F['t']])
            hits['hk_min'] += bool(c2[F['hk']] in (1.05, 1.00005))
            for simi in (0, 1):
                set_logical(clib, flib, 'i_xbl:simi', simi)
                clib.blmid(c_int(ityp)); flib.__m_xblsys_MOD_blmid(byref(c_int(ityp)))
                assert globals_diff(clib, flib, BL) == [], ('blmid', case, ityp, simi)
        if case % 25 == 0:
            assert globals_diff(clib, flib) == [], ('whole state', case)
    assert min(hits[k] for k in ('us_clamp', 'de_clamp', 'hk_min')) > 5, hits


# ------------------------------------------------------------------ part 4: difference equations

def build_station(clib, flib, args, ityp):
    blprv_call(clib, flib, args)
    clib.blkin(); flib.__m_xblsys_MOD_blkin()
    clib.blvar(c_int(ityp)); flib.__m_xblsys_MOD_blvar(byref(c_int(ityp)))


def com1_from_com2(clib, flib):
    for lib, fortran in ((clib, False), (flib, True)):
        global_array(lib, 'i_xbl:com1', np.float64, (73,), fortran)[:] = global_array(
            lib, 'i_xbl:com2', np.float64, (73,), fortran)


def set_scalar(clib, flib, name, value, dtype=np.float64):
    for lib, fortran in ((clib, False), (flib, True)):
        global_array(lib, name, dtype, None, fortran)[:] = value


def test_bldif_matches_fortran(clib, flib):
    sync_globals(clib, flib)
    clamped = unclamped = 0
    for case in range(200):
        rng = np.random.default_rng(9000 + case)
        set_bl_parameters(clib, flib, rng)
        set_scalar(clib, flib, 'i_xbl:amcrit', float(rng.uniform(4.0, 12.0)))
        set_scalar(clib, flib, 'i_xbl:idampv', int(rng.integers(0, 2)), np.int32)
        set_scalar(clib, flib, 'i_xbl:bule', float(rng.uniform(0.0, 1.0)))
        for ityp in (0, 1, 2, 3):
            styp = max(ityp, 1)
            a1 = station(rng, styp)
            build_station(clib, flib, a1, styp)
            com1_from_com2(clib, flib)
            a2 = list(station(rng, styp))
            a2[0] = a1[0] + float(rng.uniform(1e-4, 0.2))          # x2 > x1
            build_station(clib, flib, a2, styp)
            set_scalar(clib, flib, 'i_xbl:simi', int(ityp == 0), np.int32)
            clib.blmid(c_int(styp)); flib.__m_xblsys_MOD_blmid(byref(c_int(styp)))
            assert globals_diff(clib, flib, BL) == [], ('setup', case, ityp)
            clib.bldif(c_int(ityp)); flib.__m_xblsys_MOD_bldif(byref(c_int(ityp)))
            assert globals_diff(clib, flib, BL) == [], ('bldif', case, ityp)
            c1 = global_array(clib, 'i_xbl:com1', np.float64, (73,), False)
            c2 = global_array(clib, 'i_xbl:com2', np.float64, (73,), False)
            hl2 = np.log(abs((c2[F['hk']] - 1.0) / (c1[F['hk']] - 1.0))) ** 2
            clamped += hl2 > 15.0
            unclamped += hl2 <= 15.0
        if case % 40 == 0:
            assert globals_diff(clib, flib) == [], ('whole state', case)
    assert clamped > 20 and unclamped > 200, (clamped, unclamped)


def test_tesys_matches_fortran(clib, flib):
    sync_globals(clib, flib)
    for case in range(200):
        rng = np.random.default_rng(12000 + case)
        set_bl_parameters(clib, flib, rng)
        build_station(clib, flib, station(rng, 3), 3)
        cte, tte, dte = float(rng.uniform(0, 0.3)), 10.0 ** rng.uniform(-5, -2), 10.0 ** rng.uniform(-5, -2)
        clib.tesys(c_double(cte), c_double(tte), c_double(dte))
        flib.__m_xblsys_MOD_tesys(byref(c_double(cte)), byref(c_double(tte)), byref(c_double(dte)))
        assert globals_diff(clib, flib, BL) == [], case


# ------------------------------------------------------------------ part 5: transition and assembly

def laminar_interval(clib, flib, rng):
    """Stations 1 and 2 of a laminar interval with N1 below, near or above Ncrit."""
    amcrit = float(rng.uniform(4.0, 11.0))
    set_scalar(clib, flib, 'i_xbl:amcrit', amcrit)
    set_scalar(clib, flib, 'i_xbl:idampv', int(rng.integers(0, 2)), np.int32)
    a1 = list(station(rng, 1))
    a1[3] = 10.0 ** rng.uniform(-4, -2.5)                       # theta: Rtheta large enough to amplify
    a1[4] = float(rng.uniform(2.2, 3.5)) * a1[3]
    a1[1] = amcrit - float(rng.choice([rng.uniform(0, 0.05), rng.uniform(0.05, 2), rng.uniform(2, 9)]))
    build_station(clib, flib, a1, 1)
    com1_from_com2(clib, flib)
    a2 = list(station(rng, 1))
    a2[0] = a1[0] + float(rng.uniform(1e-3, 0.3))
    a2[3] = a1[3] * float(rng.uniform(1.0, 1.5))
    a2[4] = float(rng.uniform(2.2, 3.8)) * a2[3]
    a2[6] = a1[6] * float(rng.uniform(0.9, 1.1))
    build_station(clib, flib, a2, 1)
    # forced transition: none / inside the interval / at the station 2
    xi = float(rng.choice([a2[0] + 1.0, rng.uniform(a1[0], a2[0]), a2[0]]))
    set_scalar(clib, flib, 'i_xbl:xiforc', xi)
    return a1, a2


def call_trchek2(clib, flib):
    clib.trchek2.restype = c_int
    flib.__m_xblsys_MOD_trchek2.restype = c_int
    return clib.trchek2(), flib.__m_xblsys_MOD_trchek2()


def test_trchek2_and_trdif_match_fortran(clib, flib):
    sync_globals(clib, flib)
    outcomes = {'none': 0, 'free': 0, 'forced': 0}
    for case in range(400):
        rng = np.random.default_rng(15000 + case)
        set_bl_parameters(clib, flib, rng)
        laminar_interval(clib, flib, rng)
        assert globals_diff(clib, flib, BL) == [], ('setup', case)
        rc, rf = call_trchek2(clib, flib)
        assert rc == rf == 1
        assert globals_diff(clib, flib, BL) == [], ('trchek2', case)
        flags = {f: global_array(clib, 'i_xbl:' + f, np.int32, None, False)[0] for f in ('tran', 'trfree', 'trforc')}
        outcomes['forced' if flags['trforc'] else 'free' if flags['trfree'] else 'none'] += 1
        if flags['tran']:
            clib.trdif(); flib.__m_xblsys_MOD_trdif()
            assert globals_diff(clib, flib, BL) == [], ('trdif', case)
        if case % 50 == 0:
            assert globals_diff(clib, flib) == [], ('whole state', case)
    assert min(outcomes.values()) > 30, outcomes


def test_trchek_wrapper(clib, flib):
    sync_globals(clib, flib)
    rng = np.random.default_rng(1)
    set_bl_parameters(clib, flib, rng)
    laminar_interval(clib, flib, rng)
    clib.trchek.restype = c_int
    flib.__m_xblsys_MOD_trchek.restype = c_int
    assert clib.trchek() == flib.__m_xblsys_MOD_trchek() == 1
    assert globals_diff(clib, flib, BL) == []


def test_trchek2_nan_matches_fortran(clib, flib):
    """NaN station data: 30 failed Newton steps, NaN check, return .false. (ABORT_ON_NAN)."""
    sync_globals(clib, flib)
    rng = np.random.default_rng(2)
    set_bl_parameters(clib, flib, rng)
    laminar_interval(clib, flib, rng)
    for lib, fortran in ((clib, False), (flib, True)):
        global_array(lib, 'i_xbl:com2', np.float64, (73,), fortran)[F['t']] = np.nan
    rc, rf = call_trchek2(clib, flib)
    assert rc == rf == 0
    assert globals_diff(clib, flib, BL) == []


def test_blsys_matches_fortran(clib, flib):
    sync_globals(clib, flib)
    combos = 0
    for case in range(160):
        rng = np.random.default_rng(20000 + case)
        set_bl_parameters(clib, flib, rng)
        set_scalar(clib, flib, 'i_xbl:bule', float(rng.uniform(0.0, 1.0)))
        mode = case % 5            # 0 laminar, 1 turbulent, 2 wake, 3 transition, 4 similarity
        if mode == 3:
            laminar_interval(clib, flib, rng)
            call_trchek2(clib, flib)
            if not global_array(clib, 'i_xbl:tran', np.int32, None, False)[0]:
                continue
        else:
            styp = {0: 1, 1: 2, 2: 3, 4: 1}[mode]
            a1 = station(rng, styp)
            build_station(clib, flib, a1, styp)
            com1_from_com2(clib, flib)
            a2 = list(station(rng, styp))
            a2[0] = a1[0] + float(rng.uniform(1e-4, 0.2))
            build_station(clib, flib, a2, styp)
            set_scalar(clib, flib, 'i_xbl:amcrit', 9.0)
        flags = {'wake': mode == 2, 'turb': mode in (1, 2), 'tran': mode == 3, 'simi': mode == 4}
        for f, v in flags.items():
            set_scalar(clib, flib, 'i_xbl:' + f, int(v), np.int32)
        assert globals_diff(clib, flib, BL) == [], ('setup', case)
        clib.blsys(); flib.__m_xblsys_MOD_blsys()
        assert globals_diff(clib, flib, BL) == [], ('blsys', case, mode)
        combos += 1
    assert combos > 130
