# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""Phase 1: libxfoil_c builds, loads through ctypes, and passes doubles by pointer."""
from ctypes import byref, c_double

import xfoil_c
from xfoil_c import _lib


def test_library_loads():
    assert _lib.lib_path.endswith(('libxfoil_c.so', 'libxfoil_c.dylib', 'xfoil_c.dll', 'libxfoil_c.dll'))


def test_version_matches_package():
    assert xfoil_c.lib_version() == xfoil_c.__version__


def test_set_print_roundtrip(clib):
    from ctypes import c_bool
    for value in (True, False):
        clib.set_print(byref(c_bool(value)))
        assert clib.get_print() is value
    # clib fixture default: printing off (keeps C output out of the test log)


def test_context_function_list_matches_header():
    import os, re
    h = open(os.path.join(os.path.dirname(__file__), '..', 'src', 'xfoil_c.h')).read()
    in_header = set(re.findall(r'XFOIL_API\s+[\w\s\*]+?\b(\w+)\s*\(\s*xfoil_t\s*\*\s*xf\s*[,)]', h)) - {'xfoil_free'}
    assert in_header == set(_lib.CONTEXT_FUNCTIONS)


def test_instances_are_independent():
    """Each XFoil() owns its state (clean-up C5b): another instance working in between changes
    nothing, bit for bit (the viscous solution of an instance carries over between its own calls)."""
    from xfoil_c import XFoil
    from xfoil_c.model import Airfoil
    from test_api_analysis import naca0012

    def run(interleave):
        a = XFoil()
        a.print = False
        a.airfoil = Airfoil(*naca0012())
        a.Re, a.max_iter = 1e6, 100
        first = a.a(10)
        if interleave:
            b = XFoil()
            b.print = False
            b.naca('4415')
            b.Re, b.M = 3e5, 0.4
            b.a(4)
            assert (a.Re, a.M, b.Re, b.M) == (1e6, 0.0, 3e5, 0.4)
            del b
        return first, a.a(10), a.cl(0.5), a.Re, a.M

    alone, interleaved = run(False), run(True)
    assert repr(alone) == repr(interleaved)
    assert [round(v, 4) for v in alone[0][:4]] == [1.0809, 0.0150, 0.0053, -5.5077]


def test_instances_in_parallel_threads():
    """Instances in different threads run truly in parallel (ctypes releases the GIL during the
    C calls) and give bit-identical results to a sequential run (clean-up C5d)."""
    import threading
    from xfoil_c import XFoil

    def work(spec, re, out, key):
        xf = XFoil()
        xf.print = False
        xf.naca(spec)
        xf.Re, xf.max_iter = re, 100
        out[key] = repr(xf.aseq(-4, 10, 1.0)) + repr(xf.cseq(0.1, 1.0, 0.1)) + repr(xf.get_cp_distribution())

    jobs = [('2412', 1e6), ('0012', 5e5), ('4415', 3e5), ('23012', 2e6)] * 2
    sequential, parallel = {}, {}
    for k, (spec, re) in enumerate(jobs):
        work(spec, re, sequential, k)
    threads = [threading.Thread(target=work, args=(spec, re, parallel, k)) for k, (spec, re) in enumerate(jobs)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert parallel == sequential
