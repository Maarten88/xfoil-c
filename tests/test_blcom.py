# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""i_xbl COM1/COM2 aliases: after Fortran's preptrs, every pointer (X1, CF2_D2, ...) refers to
the slot of com1/com2 that the generated blcom_t field occupies."""
import json
import re
from ctypes import c_char, c_double, c_void_p, POINTER, cast
from pathlib import Path

import numpy as np

from conftest import require_xfoil_python  # noqa: E402

SRC = Path(require_xfoil_python('src'))
HDR = Path(__file__).resolve().parents[1] / 'src' / 'ixbl.h'


def test_pointer_aliases_match_blcom_fields(flib):
    flib.__m_xbl_MOD_preptrs()
    com = {k: np.frombuffer((c_char * (73 * 8)).in_dll(flib, '__i_xbl_MOD_' + k), np.float64) for k in ('com1', 'com2')}
    base = {k: com[k].ctypes.data for k in com}
    # field order in the generated header: "double <field>; /* NAME1 / NAME2 = COM(k) */"
    rows = re.findall(r'double (\w+);\s*/\* (\w+) / (\w+) = COM\((\d+)\) \*/', HDR.read_text())
    assert len(rows) == 73
    for field, n1, n2, k in rows:
        for com_name, var in (('com1', n1), ('com2', n2)):
            # gfortran pointer to a scalar: the symbol holds the target address
            addr = c_void_p.in_dll(flib, '__i_xbl_MOD_' + var.lower()).value
            assert addr == base[com_name] + 8 * (int(k) - 1), (field, var)
