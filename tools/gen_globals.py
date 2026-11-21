"""Generate the C state declarations from the Fortran global-state modules.

Usage: python -I tools/gen_globals.py <fortran_src_dir> <c_src_dir>

Writes to <c_src_dir>:
    ixfoil.h, iblpar.h, ixbl.h  parameters of i_xfoil (+ i_pindex), i_blpar, i_xbl; ixbl.h also
                                defines blcom_t (COMMON-style pointer aliases into com1/com2)
    xstate.h/.c                 typedef struct xfoil xfoil_t: all module variables as fields
                                (clean-up C5a), xfoil_state_init() (Fortran DATA values, index
                                variables -1) and xfoil_default() (the library's default context)
    xglobals.h/.c               lookup table over the fields for the tests
                                (`xfoil_global("name")` for i_xfoil, `xfoil_global("i_xbl:name")`
                                otherwise; addresses in the default context)
i_circle (inverse design) is not generated, and the variables in dropped_globals.txt (never read
by the library, clean-up C3) are left out.

Layout (user decision, Phase 4): Fortran column-major. Fortran X(d1, ..., dn) becomes
C x[dn]...[d1], so X(i,j) is xf->x[j-1][i-1] and memory is byte-identical to gfortran's.
real -> double, integer -> int, logical -> int (gfortran default logical is 4 bytes),
character(L) -> char[L] (blank-padded, no NUL). A dimension lo:hi becomes [hi - lo + 1], so
C index = Fortran index - lo. Field names are the lower-case Fortran names (a C keyword would
get a trailing underscore). Variables listed in index_globals.txt hold 0-based index values
(-1 = none) and start as -1.

i_xbl's `real, pointer` variables are aliases set by preptrs (m_xbl.f90), e.g.
X1 => COM1(1), X2 => COM2(1). Slot k means the same quantity in both arrays, so com1
and com2 are a union `blcom_t` of a double array a[NCOM] and a struct of named fields:
Fortran X1 is com1.x, CF2_D2 is com2.cf_d, COM1(k) is com1.a[k-1].
"""
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

src, out = Path(sys.argv[1]), Path(sys.argv[2])
HERE = Path(__file__).resolve().parent

MODULES = [('i_xfoil', 'ixfoil', ['i_pindex']), ('i_blpar', 'iblpar', []), ('i_xbl', 'ixbl', [])]
INCLUDES = {'ixfoil': [], 'iblpar': [], 'ixbl': ['iblpar.h']}


def joined_lines(path):
    """Statements with continuation lines joined; comments returned separately."""
    stmts, comments, buf = [], [], ''
    for raw in path.read_text(errors='replace').splitlines():
        code, _, comment = raw.partition('!')
        if comment.strip():
            comments.append(comment)
        code = code.rstrip()
        if not code.strip():
            continue
        line = (buf + ' ' + code.strip().lstrip('&')) if buf else code.strip()
        if line.endswith('&'):
            buf = line[:-1]
            continue
        buf = ''
        stmts.append(line)
    return stmts, comments


def split_top(body):
    items, depth, cur = [], 0, ''
    for ch in body:
        depth += ch == '('
        depth -= ch == ')'
        if ch == ',' and depth == 0:
            items.append(cur.strip())
            cur = ''
        else:
            cur += ch
    items.append(cur.strip())
    return items


DECL = re.compile(r'^(real|integer|logical|character|complex)\s*(\(\s*(\d+)\s*\))?\s*'
                  r'((?:,\s*[a-z]+(?:\([^)]*\))?)*)\s*::\s*(.*)$', re.I)


def parse(path, module):
    stmts, comments = joined_lines(path)
    params, variables, pointers = [], [], []
    for st in stmts:
        m = DECL.match(st)
        if not m:
            continue
        typ, clen, attrs, body = m.group(1).lower(), m.group(3), m.group(4).lower(), m.group(5)
        dim_attr = re.search(r'dimension\s*\(([^)]*)\)', attrs)
        for it in split_top(body):
            im = re.match(r'([a-z0-9_]+)\s*(\(([^)]*)\))?\s*(=\s*(.*))?$', it, re.I)
            if not im:
                raise SystemExit('cannot parse declaration item %r in %s' % (it, path))
            name, dims, init = im.group(1), im.group(3) or (dim_attr.group(1) if dim_attr else None), im.group(5)
            if 'parameter' in attrs:
                params.append((name.upper(), init.strip(), typ))
            elif 'pointer' in attrs:
                pointers.append(name.lower())
            else:
                variables.append(dict(module=module, name=name.lower(), fname=name.upper(), type=typ,
                                      clen=int(clen) if clen else None,
                                      dims=[d.strip() for d in split_top(dims)] if dims else [],
                                      init=init.strip() if init else None))
    docs = {}
    for c in comments:
        dm = re.match(r'\s+([A-Z][A-Z0-9_]*(?:\s*,\s*[A-Z][A-Z0-9_]*)*)\s{2,}(\S.*)$', c)
        if dm:
            for nm in re.split(r'\s*,\s*', dm.group(1)):
                docs.setdefault(nm.lower(), dm.group(2).strip())
    return params, variables, pointers, docs


# ---- parse all modules
mods = {}
all_params = []
for module, base, extra in MODULES:
    params = []
    for e in extra:
        params += parse(src / (e + '.f90'), e)[0]
    p, v, ptrs, docs = parse(src / (module + '.f90'), module)
    params += p
    mods[module] = dict(base=base, params=params, variables=v, pointers=ptrs, docs=docs)
    all_params += params
pnames = {n for n, _, _ in all_params}
variables = [v for m in mods.values() for v in m['variables']]


def c_expr(e):
    """Fortran integer constant expression -> C (parameter names upper-cased)."""
    return re.sub(r'[A-Za-z_][A-Za-z0-9_]*', lambda m: m.group(0).upper() if m.group(0).upper() in pnames
                  else m.group(0), e)


# ---- i_xbl COMMON aliases (preptrs in m_xbl.f90)
xbl_src = (src / 'm_xbl.f90').read_text(errors='replace')
pre = re.search(r'subroutine\s+preptrs(.*?)end\s+subroutine\s+preptrs', xbl_src, re.I | re.S).group(1)
assoc = re.findall(r'^\s*([A-Za-z0-9_]+)\s*=>\s*(COM[12])\s*\(\s*(\d+)\s*\)', pre, re.M)
slots = {('com1', int(k)): n.lower() for n, c, k in assoc if c == 'COM1'}
slots.update({('com2', int(k)): n.lower() for n, c, k in assoc if c == 'COM2'})
NCOM = int(dict((n, e) for n, e, _ in all_params)['NCOM'])
fields = []
for k in range(1, NCOM + 1):
    a, b = slots[('com1', k)], slots[('com2', k)]
    if len(a) != len(b) or any(x != y and not (x == '1' and y == '2') for x, y in zip(a, b)):
        raise SystemExit('COM slot %d: %s / %s is not a name1/name2 pair' % (k, a, b))
    fields.append((''.join(x for x, y in zip(a, b) if x == y), a, b))
if sorted(set(slots.values())) != sorted(mods['i_xbl']['pointers']):
    raise SystemExit('i_xbl pointers and preptrs associations differ')
if len({f for f, _, _ in fields}) != NCOM:
    raise SystemExit('COM field names not unique')
BLCOM = {'com1', 'com2', 'c1sav', 'c2sav'}   # c1sav/c2sav: saved copies of com1/com2 (trdif, trchek2)

# ---- globals holding index values: 0-based in C, -1 = none (see index_globals.txt)
INDEX = {l.split('#', 1)[0].strip() for l in (HERE / 'index_globals.txt').read_text().splitlines()} - {''}
INDEX = {i if ':' in i else 'i_xfoil:' + i for i in INDEX}
_vars = {v['module'] + ':' + v['name']: v for v in variables}
for key in sorted(INDEX):
    if key not in _vars or _vars[key]['type'] != 'integer':
        raise SystemExit('index_globals.txt: %s is not an integer global' % key)

# ---- field names: the Fortran name (lower case); a C keyword would get a trailing underscore
C_KEYWORDS = set('''auto break case char const continue default do double else enum extern float for goto if
inline int long register restrict return short signed sizeof static struct switch typedef union unsigned void
volatile while bool true false'''.split())
for v in variables:
    v['cname'] = v['name'] + '_' if v['name'] in C_KEYWORDS else v['name']


# ---- globals the library never reads (clean-up C3): not generated. Done after the renames so
# the remaining names do not change.
DROP = {l.split('#', 1)[0].strip() for l in (HERE / 'dropped_globals.txt').read_text().splitlines()} - {''}
DROP = {d if ':' in d else 'i_xfoil:' + d for d in DROP}
for key in sorted(DROP):
    if key not in _vars:
        raise SystemExit('dropped_globals.txt: unknown global %s' % key)
    if key in INDEX:
        raise SystemExit('dropped_globals.txt: %s is also in index_globals.txt' % key)
for m in mods.values():
    m['variables'] = [v for v in m['variables'] if v['module'] + ':' + v['name'] not in DROP]
variables = [v for m in mods.values() for v in m['variables']]


def ctype(v):
    if v['module'] == 'i_xbl' and v['name'] in BLCOM:
        return 'blcom_t'
    return {'real': 'double', 'integer': 'int', 'logical': 'int', 'character': 'char',
            'complex': 'double _Complex'}[v['type']]


def c_dim(d):
    """Fortran extent 'n' or 'lo:hi' -> C array size."""
    if ':' in d:
        lo, hi = (x.strip() for x in d.split(':'))
        return '%s + 1' % c_expr(hi) if lo == '0' else '%s - %s + 1' % (c_expr(hi), c_expr(lo))
    return c_expr(d)


def cdims(v):
    if v['module'] == 'i_xbl' and v['name'] in BLCOM:
        return ''
    dims = ''.join('[%s]' % c_dim(d) for d in reversed(v['dims']))
    if v['type'] == 'character':
        dims += '[%d]' % v['clen']
    return dims


def cinit(v):
    if v['init'] is None:
        return ''
    i = v['init'].lower()
    if v['type'] == 'logical':
        return ' = %d' % (1 if i == '.true.' else 0)
    return ' = %s' % i


def fdecl(v):
    t = v['type'] + ('(%d)' % v['clen'] if v['clen'] else '')
    d = '(%s)' % ', '.join(v['dims']) if v['dims'] else ''
    return '%s %s%s' % (t, v['fname'], d)


def banner(module):
    return '''/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  GENERATED by tools/gen_globals.py from %s.f90. Do not edit by hand;
 *  regenerate instead.
 *
 *  Arrays keep the Fortran (column-major) memory layout: Fortran X(i, j) is C
 *  x[j-1][i-1] (dimensions reversed, 0-based). Types: real -> double,
 *  integer -> int, logical -> int (4 bytes, 1 = .true.), character(L) -> char[L]
 *  (blank-padded, no NUL).
 */
''' % module


for module, base, _ in MODULES:
    m = mods[module]
    guard = 'XFOIL_C_%s_H' % base.upper()
    h = [banner(module), '#ifndef %s\n#define %s\n\n' % (guard, guard)]
    h += ['#include "%s"\n' % inc for inc in INCLUDES[base]]
    if INCLUDES[base]:
        h.append('\n')
    int_params = [(n, e, t) for n, e, t in m['params'] if t in ('integer', 'logical')]
    if int_params:
        h.append('/* ---- parameters (%s) ---- */\nenum {\n' % ', '.join([module] + [e for mm, b, ex in MODULES
                                                                                       if mm == module for e in ex]))
        for n, e, t in int_params:
            if t == 'integer':
                h.append('    %s = %s,\n' % (n, c_expr(e)))
            else:
                h.append('    %s = %d,   /* logical %s */\n' % (n, 1 if e.lower() == '.true.' else 0, e))
        h.append('};\n')
    for n, e, t in m['params']:
        if t == 'real':
            h.append('static const double %s = %s;\n' % (n, c_expr(e)))
        elif t not in ('integer', 'logical'):
            raise SystemExit('unsupported parameter type %s for %s' % (t, n))
    h.append('\n')
    if module == 'i_xbl':
        h.append('/* One BL station of the COMMON blocks: Fortran COM1/COM2(NCOM) with the pointer\n'
                 ' * aliases set by preptrs (m_xbl.f90). Fortran X1 is com1.x, X2 is com2.x,\n'
                 ' * COM1(k) is com1.a[k-1]. */\ntypedef union {\n    double a[NCOM];\n    struct {\n')
        for k, (f, a, b) in enumerate(fields, 1):
            h.append('        double %-8s /* %s / %s = COM(%d) */\n' % (f + ';', a.upper(), b.upper(), k))
        h.append('    };\n} blcom_t;\n_Static_assert(sizeof(blcom_t) == NCOM * sizeof(double), "blcom_t layout");\n\n')
    h.append('/* The variables of %s are fields of xfoil_t (xstate.h). */\n' % module)
    h.append('\n#endif /* %s */\n' % guard)
    (out / (base + '.h')).write_text(''.join(h))
    stale = out / (base + '.c')
    if stale.exists():
        stale.unlink()

# ---- the state struct (xstate.h/.c): all module variables, one instance per XFoil context
sh = ['''/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  GENERATED by tools/gen_globals.py from the Fortran state modules. Do not edit by hand.
 *
 *  xfoil_t holds the whole state of one XFoil instance: the variables of the Fortran
 *  modules i_xfoil, i_blpar and i_xbl (minus tools/dropped_globals.txt) as fields with the
 *  Fortran names. Routines reach it through their first parameter `xfoil_t *xf`.
 *  Arrays keep the Fortran (column-major) memory layout: Fortran X(i, j) is C
 *  xf->x[j-1][i-1]. Types: real -> double, integer -> int, logical -> int (1 = .true.),
 *  character(L) -> char[L] (blank-padded, no NUL).
 */
#ifndef XFOIL_C_XSTATE_H
#define XFOIL_C_XSTATE_H

#include "xfoil_c.h"
#include "ixfoil.h"
#include "iblpar.h"
#include "ixbl.h"
#include "xscratch.h"

typedef struct xfoil {
''']
for module, base, _ in MODULES:
    m = mods[module]
    sh.append('    /* ---- %s ---- */\n' % module)
    for v in m['variables']:
        doc = m['docs'].get(v['name'])
        mark = ' [0-based index values; -1 = none]' if module + ':' + v['name'] in INDEX else ''
        sh.append('    /* %s%s%s */\n    %s %s%s;\n' % (fdecl(v), (': ' + doc.replace('*/', '* /')) if doc else '', mark,
                                                   ctype(v), v['cname'], cdims(v)))
sh.append('''    /* ---- work arrays of single routines (not Fortran module variables; xscratch.h) ---- */
    xfoil_scratch_t scr;
} xfoil_t;

/* Set a ZEROED context (calloc) to the Fortran initial state (DATA values, index variables -1). */
void xfoil_state_init(xfoil_t *xf);

/* The library's default context (allocated on first use). */
XFOIL_API xfoil_t *xfoil_default(void);

/* Print flag of the running API call (the context's show_output, set on entry). */
extern _Thread_local int xfoil_print;

#endif /* XFOIL_C_XSTATE_H */
''')
(out / 'xstate.h').write_text(''.join(sh))
sc = ['''/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  GENERATED by tools/gen_globals.py. Do not edit by hand.
 */
#include "xstate.h"

#include <stdlib.h>
#include <string.h>
#include <threads.h>

_Thread_local int xfoil_print = 1;

void xfoil_state_init(xfoil_t *xf)
{
''']
for v in variables:
    if v['init'] is not None:
        if v['dims'] or v['type'] == 'character':
            raise SystemExit('initialised array/string not supported: %s' % v['name'])
        sc.append('    xf->%s%s;\n' % (v['cname'], cinit(v)))
for k in sorted(INDEX):
    sc.append('    memset(&xf->%s, 0xff, sizeof xf->%s);   /* all -1 */\n' % (_vars[k]['cname'], _vars[k]['cname']))
sc.append('''}

static xfoil_t *def = NULL;
static once_flag def_once = ONCE_FLAG_INIT;

static void def_create(void)
{
    def = calloc(1, sizeof *def);
    if (!def) abort();
    xfoil_state_init(def);
}

xfoil_t *xfoil_default(void)
{
    call_once(&def_once, def_create);     /* thread-safe creation on first use */
    return def;
}
''')
(out / 'xstate.c').write_text(''.join(sc))

# ---- lookup table over all modules + index-global initialisation
xh = '''/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  GENERATED by tools/gen_globals.py. Lookup of the Fortran-module globals, for
 *  the tests (which compare them with the Fortran symbols __<module>_MOD_<name>).
 */
#ifndef XFOIL_C_XGLOBALS_H
#define XFOIL_C_XGLOBALS_H

#include <stddef.h>

#include "xfoil_c.h"

/* Address and size in bytes of a global, or NULL if unknown. `name` is the lower-case
 * Fortran name, "module:name" (e.g. "i_xbl:amcrit") or just "name" for i_xfoil. */
XFOIL_API void *xfoil_global(const char *name, size_t *size);

/* Number of globals, the "module:name" of entry k, and its Fortran type:
 * 'r' real (double), 'i' integer, 'l' logical, 'c' character, 'x' complex (2 doubles). */
XFOIL_API int xfoil_global_count(void);
XFOIL_API const char *xfoil_global_name(int k);
XFOIL_API char xfoil_global_type(int k);

#endif /* XFOIL_C_XGLOBALS_H */
'''
(out / 'xglobals.h').write_text(xh)
xc = [xh.split('#ifndef')[0].replace('Lookup of the Fortran-module globals', 'Lookup table of the Fortran-module globals'),
      '#include "xglobals.h"\n\n#include <string.h>\n\n']
xc.append('#include "xstate.h"\n')
xc.append('\nstatic const struct {\n    const char *name;     /* module:name */\n    size_t offset;        /* in xfoil_t */\n    size_t size;\n'
          '    char type;            /* r real, i integer, l logical, c character, x complex */\n} TABLE[] = {\n')
for v in variables:
    code = 'x' if v['type'] == 'complex' else v['type'][0]
    xc.append('    {"%s:%s", offsetof(xfoil_t, %s), sizeof ((xfoil_t *)0)->%s, \'%s\'},\n' % (v['module'], v['name'], v['cname'], v['cname'], code))
xc.append('''};

enum { NGLOBALS = (int)(sizeof TABLE / sizeof TABLE[0]) };

void *xfoil_global(const char *name, size_t *size)
{
    char key[96];
    if (!strchr(name, ':')) {               /* bare name: i_xfoil */
        strcpy(key, "i_xfoil:");
        strncat(key, name, sizeof key - strlen(key) - 1);
        name = key;
    }
    for (int k = 0; k < NGLOBALS; k++) {
        if (strcmp(TABLE[k].name, name) == 0) {
            if (size) *size = TABLE[k].size;
            return (char *)xfoil_default() + TABLE[k].offset;
        }
    }
    return NULL;
}

int xfoil_global_count(void)
{
    return NGLOBALS;
}

const char *xfoil_global_name(int k)
{
    return (k >= 0 && k < NGLOBALS) ? TABLE[k].name : NULL;
}

char xfoil_global_type(int k)
{
    return (k >= 0 && k < NGLOBALS) ? TABLE[k].type : 0;
}
''')
(out / 'xglobals.c').write_text(''.join(xc))

print('modules: %s' % ', '.join('%s (%d vars%s)' % (m, len(mods[m]['variables']),
                                                    ', %d pointer aliases' % len(mods[m]['pointers'])
                                                    if mods[m]['pointers'] else '') for m, _, _ in MODULES))
print('parameters: %d' % len(all_params))
print('index globals: %s' % ', '.join(sorted(INDEX)))
print('dropped: %d' % len(DROP))
