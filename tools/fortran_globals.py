"""Parse a Fortran include module's declarations and report which of its variables
each routine of another module references.

Usage: python -I tools/fortran_globals.py <src_dir> <i_module.f90> <m_module.f90> [--reachable-only] [--json out]
"""
import json
import re
import subprocess
import sys
from pathlib import Path

src = Path(sys.argv[1])
inc_file, mod_file = sys.argv[2], sys.argv[3]
reach_only = '--reachable-only' in sys.argv

# ---- parameters and variable declarations of the include module (joined continuations)
lines, buf = [], ''
for raw in (src / inc_file).read_text(errors='replace').splitlines():
    line = raw.split('!', 1)[0].rstrip()
    if not line.strip():
        continue
    if buf:
        line = buf + ' ' + line.strip().lstrip('&')
    if line.endswith('&'):
        buf = line[:-1]
        continue
    buf = ''
    lines.append(line.strip())

params, decls = {}, {}
DECL = re.compile(r'^(real|integer|logical|character|complex)\s*(\([^)]*\))?\s*((?:,\s*[a-z]+(?:\([^)]*\))?)*)\s*::\s*(.*)$', re.I)
for line in lines:
    m = DECL.match(line)
    if not m:
        continue
    typ = (m.group(1) + (m.group(2) or '')).lower()
    attrs = m.group(3).lower()
    body = m.group(4)
    # split on commas not inside parentheses
    items, depth, cur = [], 0, ''
    for ch in body:
        depth += ch == '('
        depth -= ch == ')'
        if ch == ',' and depth == 0:
            items.append(cur.strip()); cur = ''
        else:
            cur += ch
    items.append(cur.strip())
    for it in items:
        name_m = re.match(r'([a-z0-9_]+)\s*(\(([^)]*)\))?\s*(=\s*(.*))?$', it, re.I)
        if not name_m:
            continue
        name = name_m.group(1).lower()
        dims = name_m.group(3)
        init = name_m.group(5)
        if 'parameter' in attrs:
            params[name] = init
        else:
            dm = re.search(r'dimension\s*\(([^)]*)\)', attrs)
            decls[name] = dict(type=typ, dims=(dims or (dm.group(1) if dm else None)), init=init)

# evaluate parameter values
env = {}
for _ in range(5):
    for k, v in params.items():
        try:
            env[k] = int(eval(v.lower().replace('/', '//'), {}, env))
        except Exception:
            pass

def size_of(d):
    if not d['dims']:
        return 1
    n = 1
    for part in d['dims'].lower().split(','):
        lo, _, hi = part.partition(':')
        lo, hi = (lo, hi) if hi else ('1', lo)
        n *= int(eval(hi.replace('/', '//'), {}, env)) - int(eval(lo.replace('/', '//'), {}, env)) + 1
    return n

for d in decls.values():
    try:
        d['size'] = size_of(d)
    except Exception:
        d['size'] = None

# ---- routines of the target module and their references
# call graph (to a temp json) for reachability
import tempfile, os
tmp = tempfile.mktemp(suffix='.json')
subprocess.run([sys.executable, '-I', str(Path(__file__).with_name('fortran_callgraph.py')), str(src), '--json', tmp],
               capture_output=True, text=True, check=True)
cg = json.load(open(tmp)); os.remove(tmp)
reach = set(cg['reachable'])
mod_lines = (src / mod_file).read_text(errors='replace').splitlines()
usage = {}
for name, r in cg['routines'].items():
    if r['file'] != mod_file or (reach_only and name not in reach):
        continue
    text = '\n'.join(l.split('!', 1)[0].lower() for l in mod_lines[r['start'] - 1:r['end']])
    used = sorted(v for v in decls if re.search(r'(?<![a-z0-9_%])' + v + r'(?![a-z0-9_])', text))
    usage[name] = used

allused = sorted(set(v for u in usage.values() for v in u))
print('parameters:', ', '.join('%s=%s' % kv for kv in env.items()))
print('%s declares %d variables (%d doubles total)' % (inc_file, len(decls),
      sum(d['size'] or 0 for d in decls.values() if d['type'].startswith('real'))))
print('%s routines analysed: %d; distinct %s variables used: %d' % (mod_file, len(usage), inc_file, len(allused)))
for name in sorted(usage, key=lambda n: cg['routines'][n]['start']):
    print('  %-8s (%3d) %s' % (name, len(usage[name]), ' '.join(usage[name])))
big = sorted(((decls[v]['size'] or 0, v, decls[v]['type'], decls[v]['dims']) for v in allused), reverse=True)[:12]
print('largest used:', '; '.join('%s %s(%s)=%d' % (t, v, d, s) for s, v, t, d in big))
if '--json' in sys.argv:
    json.dump({'params': env, 'decls': decls, 'usage': usage}, open(sys.argv[sys.argv.index('--json') + 1], 'w'), indent=1)
