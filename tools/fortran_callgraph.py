"""Static call graph of the Fortran library (excluding p_xfoil.f90).

Usage: python callgraph.py <src_dir> [--json out.json]
Prints the routines reachable from the Python-facing API entry points.
Function references are matched as `name(` for known function names, so the
result may over-estimate reachability slightly (never under-estimate calls).
"""
import json
import re
import sys
from pathlib import Path

src = Path(sys.argv[1])
HDR = re.compile(r'^\s*(?:(?:recursive|pure|elemental|real|integer|logical|character\S*)\s+)*'
                 r'(subroutine|function)\s+([a-z0-9_]+)', re.I)
END = re.compile(r'^\s*(?:\d+\s+)?end\s+(subroutine|function)\b', re.I)
BIND = re.compile(r"bind\s*\(\s*c\s*,\s*name\s*=\s*'([a-z0-9_]+)'", re.I)

routines = {}       # name -> dict(file, kind, start, end, body)
for f in sorted(src.glob('*.f90')):
    if f.name == 'p_xfoil.f90':
        continue
    cur = None
    lines = f.read_text(errors='replace').splitlines()
    for lineno, raw in enumerate(lines, 1):
        line = raw.split('!', 1)[0]
        m = HDR.match(line)
        if m and cur is None:
            name = m.group(2).lower()
            hdr, k = line, lineno
            while hdr.rstrip().endswith('&'):          # header continued on next line(s)
                hdr = hdr.rstrip()[:-1] + lines[k].split('!', 1)[0]
                k += 1
            b = BIND.search(hdr)
            cur = dict(name=name, file=f.name, kind=m.group(1).lower(), start=lineno,
                       cname=b.group(1).lower() if b else None, body=[])
            continue
        if cur is not None:
            if END.match(line):
                cur['end'] = lineno
                routines[cur['name']] = cur
                cur = None
            else:
                cur['body'].append(line.lower())

functions = {n for n, r in routines.items() if r['kind'] == 'function'}
CALL = re.compile(r'\bcall\s+([a-z0-9_]+)')
graph = {}
for n, r in routines.items():
    text = '\n'.join(r['body'])
    callees = set(CALL.findall(text))
    for fn in functions:
        if fn != n and re.search(r'(?<![a-z0-9_%])' + fn + r'\s*\(', text):
            callees.add(fn)
    graph[n] = sorted(c for c in callees if c in routines)

entry = sorted(n for n, r in routines.items() if r['cname'] and r['cname'] != 'xfoil')
seen, stack = set(), list(entry)
while stack:
    n = stack.pop()
    if n in seen:
        continue
    seen.add(n)
    stack.extend(graph.get(n, []))

by_file = {}
for n in sorted(routines):
    r = routines[n]
    by_file.setdefault(r['file'], []).append((n, n in seen, r['end'] - r['start'] + 1))

print('API entry points (%d): %s' % (len(entry), ' '.join(entry)))
print('Reachable: %d of %d routines' % (len(seen), len(routines)))
tot_lines = sum(r['end'] - r['start'] + 1 for r in routines.values())
reach_lines = sum(routines[n]['end'] - routines[n]['start'] + 1 for n in seen)
print('Lines in routines: reachable %d of %d' % (reach_lines, tot_lines))
print()
print('%-14s %9s %13s  unreachable routines' % ('file', 'reach/all', 'reach lines'))
for f, items in by_file.items():
    rl = sum(l for _, ok, l in items if ok)
    al = sum(l for _, _, l in items)
    unreach = [n for n, ok, _ in items if not ok]
    print('%-14s %4d/%-4d %6d/%-6d  %s' % (f, sum(ok for _, ok, _ in items), len(items), rl, al,
                                           ' '.join(unreach)))
if '--json' in sys.argv:
    out = sys.argv[sys.argv.index('--json') + 1]
    json.dump({'entry': entry, 'reachable': sorted(seen), 'graph': graph,
               'routines': {n: {k: v for k, v in r.items() if k != 'body'} for n, r in routines.items()}},
              open(out, 'w'), indent=1)
