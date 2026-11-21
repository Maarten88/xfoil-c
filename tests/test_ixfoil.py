# -*- coding: utf-8 -*-
#   This file is part of xfoil_c, a C port of XFoil.
#   Licensed under the GNU General Public License v3 or later. See LICENSE.
"""Generated module globals (src/ixfoil, iblpar, ixbl + xglobals table): same names, sizes and
initial bytes as the Fortran module variables (index globals mapped: Fortran 0 = none -> C -1).

Run in a fresh process so that no test has modified the state yet. The only Fortran symbols
without a C counterpart are i_xbl's pointer descriptors (replaced by the blcom_t union) and the
globals the library never reads (tools/dropped_globals.txt, clean-up C3; i_circle is not generated).
"""
import re
import subprocess
import sys
import textwrap
from pathlib import Path

from conftest import FORTRAN_LIB, INDEX_KEYS, TOOLS, require_xfoil_python

MODULES = ('i_xfoil', 'i_blpar', 'i_xbl')
SRC = Path(require_xfoil_python('src'))
DROPPED_FILE = Path(TOOLS) / 'dropped_globals.txt'
DROPPED = {d if ':' in d else 'i_xfoil:' + d for d in
           {l.split('#', 1)[0].strip() for l in DROPPED_FILE.read_text().splitlines()} - {''}}

SCRIPT = textwrap.dedent('''
    import subprocess, sys
    from ctypes import CDLL, c_char, c_char_p, c_int, c_size_t, c_void_p, byref
    import numpy as np
    clib = CDLL(sys.argv[1]); flib = CDLL(sys.argv[2])
    INDEX = set(sys.argv[3].split(','))
    MODULES = sys.argv[4].split(',')
    clib.xfoil_global.restype = c_void_p
    clib.xfoil_global_name.restype = c_char_p
    nm = {}
    for line in subprocess.run(['nm', '-S', '--defined-only', sys.argv[2]], capture_output=True,
                               text=True, check=True).stdout.splitlines():
        p = line.split()
        for m in MODULES:
            if len(p) == 4 and p[3].startswith('__%s_MOD_' % m):
                nm[m + ':' + p[3][len('__%s_MOD_' % m):]] = int(p[1], 16)
    names = [clib.xfoil_global_name(c_int(k)).decode() for k in range(clib.xfoil_global_count())]
    problems = []
    for name in names:
        size = c_size_t()
        addr = clib.xfoil_global(name.encode(), byref(size))
        if size.value != nm.get(name):
            problems.append('%s: C %d bytes, Fortran %s' % (name, size.value, nm.get(name)))
            continue
        module, var = name.split(':')
        cbytes = (c_char * size.value).from_address(addr).raw
        fbytes = (c_char * size.value).in_dll(flib, '__%s_MOD_%s' % (module, var)).raw
        if name in INDEX:
            f = np.frombuffer(fbytes, np.int32)
            fbytes = np.where(f >= 1, f - 1, -1).astype(np.int32).tobytes()
        if cbytes != fbytes:
            problems.append('%s: initial bytes differ' % name)
    sc = (c_char * 4).from_address(clib.xfoil_global(b'show_output', None)).raw
    print('COUNT', len(names), 'SHOW_OUTPUT', sc.hex(), 'BYTES', sum(nm.get(n, 0) for n in names))
    print('PROBLEMS', problems)
    print('FORTRAN_ONLY', ','.join(sorted(set(nm) - set(names))))
''')


def xbl_pointer_names():
    text = re.sub(r'&\s*\n\s*&?', ' ', (SRC / 'i_xbl.f90').read_text())
    names = set()
    for m in re.finditer(r'real\s*,\s*pointer\s*::(.*)', text, re.I):
        names |= {'i_xbl:' + n.strip().lower() for n in m.group(1).split(',') if n.strip()}
    return names


def test_globals_match_fortran_layout_and_initial_state(clib, flib):
    from xfoil_c import _lib
    p = subprocess.run([sys.executable, '-c', SCRIPT, _lib.lib_path, FORTRAN_LIB, ','.join(sorted(INDEX_KEYS)),
                        ','.join(MODULES)], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    first, problems, fortran_only = p.stdout.strip().splitlines()
    _, count, _, show_output, _, nbytes = first.split()
    assert problems == 'PROBLEMS []'
    assert count == str(428 + 10 + 56 - len(DROPPED)) and int(nbytes) > 25_000_000
    assert show_output == '01000000'            # logical(4) .true., little-endian
    pointers = xbl_pointer_names()
    assert len(pointers) == 146
    assert set(fortran_only.split(' ', 1)[1].split(',')) == pointers | DROPPED


def test_unknown_global(clib):
    from ctypes import c_void_p
    clib.xfoil_global.restype = c_void_p
    assert clib.xfoil_global(b'no_such_variable', None) is None
    assert clib.xfoil_global(b'i_xbl:amcrit', None) is not None
