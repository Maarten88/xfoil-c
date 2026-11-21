# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""s_xbl.f90 iblsys -> src/sxbl.c: global state identical to Fortran after the index mapping
(isys is an index global: C value = Fortran value - 1, -1 = none = Fortran 0)."""
import numpy as np

from conftest import INDEX_GLOBALS, global_array, to_c_index

IVX, ISX = 452, 2


def test_iblsys_matches_fortran(clib, flib):
    assert 'isys' in INDEX_GLOBALS
    views = {}
    for tag, lib, fortran in (('c', clib, False), ('f', flib, True)):
        views[tag] = dict(nbl=global_array(lib, 'nbl', np.int32, (ISX,), fortran),
                          isys=global_array(lib, 'isys', np.int32, (IVX, ISX), fortran),
                          nsys=global_array(lib, 'nsys', np.int32, None, fortran))
    rng = np.random.default_rng(0)
    cases = [(0, 0), (1, 1), (2, 2), (IVX, IVX), (IVX, 1)] + [tuple(rng.integers(0, IVX + 1, 2)) for _ in range(50)]
    for nbl in cases:
        views['f']['isys'][:] = 0                   # "unset" in each convention
        views['c']['isys'][:] = -1
        for tag in ('c', 'f'):
            views[tag]['nsys'][:] = -7
            views[tag]['nbl'][:] = nbl
        clib.iblsys()
        flib.__s_xbl_MOD_iblsys()
        assert np.array_equal(views['c']['isys'], to_c_index(views['f']['isys'])), nbl
        assert views['c']['nsys'].tobytes() == views['f']['nsys'].tobytes()
        assert views['c']['nbl'].tobytes() == views['f']['nbl'].tobytes()
        assert views['c']['nsys'][0] == max(nbl[0] - 1, 0) + max(nbl[1] - 1, 0)
    assert views['c']['isys'][1, 0] == 0          # first line number is 0 in C (1 in Fortran)
