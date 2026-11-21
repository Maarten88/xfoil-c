# xfoil_c

A C port of [xfoil-python](https://github.com/daniel-de-vries/xfoil-python), the Python module around Mark Drela's
XFOIL airfoil analysis code. The Fortran core is translated to C. The reference is xfoil-python with a few
changes (`docs/xfoil-python.patch`; see "Differences from upstream xfoil-python" below); its Python API is
kept, so code written for it switches with one import:

```python
from xfoil_c import XFoil          # instead of: from xfoil import XFoil
from xfoil_c.model import Airfoil
```

- **Same results.** The default (exact) build reproduces the Fortran results of the reference bit for bit (on
  the same platform; see "Tests"); this is checked routine by routine against the Fortran library (about 600
  tests).
- **Faster.** About 0.93x the Fortran run time in the exact build, about 0.67x with the optional fast build for
  your CPU (viscous analyses).
- **Independent instances.** Every `XFoil()` has its own state, and instances can run in parallel threads (the
  reference shares one Fortran state between all instances; upstream xfoil-python loads a temporary copy of
  the library per instance).
- **Errors as exceptions.** Hard errors in XFOIL (e.g. array limits) raise `xfoil_c.XFoilError` with XFOIL's
  message instead of ending the Python process.
- **Extras:** airfoil geometry properties, flap deflection and hinge moments, the boundary layer (XFOIL's DUMP).

Status: tested on Linux (x86-64, gcc 14, CMake 3.31, Python 3.12) and Windows (x86-64, MinGW-w64 gcc with the
UCRT runtime, CMake, python.org Python). Other platforms are untested.

## Install

Requirements: a C compiler (gcc), CMake >= 3.10, Python 3 with numpy.

```bash
pip install .                                   # from a checkout of this repository
```

**Windows:** build with MinGW-w64 gcc, not MSVC (the build uses gcc flags such as `-ffp-contract=off`). Use a
UCRT build of MinGW-w64 gcc (WinLibs UCRT or MSYS2 UCRT64) and put its `bin` directory and CMake on `PATH`;
setup.py then selects the "MinGW Makefiles" generator (override with `CMAKE_GENERATOR`). The DLL links libgcc
statically and needs only system DLLs (KERNEL32 and the Windows UCRT), so it loads from a python.org Python
without extra DLLs on `PATH`. Results use the Windows C runtime's math functions (`exp`, `pow`, `atan2`, ...)
and so differ from a Linux build in the last bits.

## Usage

```python
from xfoil_c import XFoil

xf = XFoil()
xf.print = False                   # XFOIL's console output off
xf.naca('2412')                    # or: xf.airfoil = Airfoil(x, y)
xf.repanel()                       # optional: XFOIL's default paneling (160 nodes)
xf.Re = 1e6
xf.max_iter = 100

cl, cd, cm, cp_min, xtr_top, xtr_bot = xf.a(4.0)        # one angle of attack (with transition points)
a, cl, cd, cm, cp_min = xf.aseq(-4, 12, 0.5)            # a sequence (end value excluded, like range())
x, y, cp = xf.get_cp_distribution()                     # Cp of the last analysed point
```

Unconverged points are returned as NaN, as in xfoil-python.

### Extensions (not in xfoil-python)

Read-only geometry of the current airfoil, computed as XFOIL's GEOPAR/TCCALC do: `max_thickness`,
`x_max_thickness`, `max_camber`, `x_max_camber`, `le_radius`, `te_angle` (degrees), `area`.

```python
xf.naca('2412')
xf.max_camber, xf.x_max_camber      # (0.0200, 0.3996)
```

Flap and hinge moment (XFOIL's GDES FLAP and OPER FMOM, without prompts):

```python
xf.set_hinge(0.75, 0.5, y_relative=True)   # hinge at x = 0.75, halfway between the surfaces
xf.a(4.0); xf.hinge_moment                # hinge moment (and hinge_fx, hinge_fy) of the last analysed point
xf.flap(5.0)                              # flap deflection 5 deg (TE down), then repanel
xf.flap(-3.0)                             # -3 deg (absolute, from the undeflected airfoil)
xf.flap(0.0)                              # back to the undeflected airfoil
```

`flap(angle)` sets the deflection (XFOIL's FLAP command adds to it) and repanels with the current paneling
settings; `flap_angle` returns it and `hinge` the hinge point. Moving the hinge keeps the deflection.
`xf.airfoil` returns the deflected airfoil. A new airfoil (`airfoil = ...`, `naca(...)`) clears hinge and deflection.

Boundary layer of the last analysed point, with the quantities of XFOIL's OPER DUMP:

```python
xf.a(6.0)
bl = xf.get_bl()                 # None unless the last point is a converged viscous one
top = bl['top']                  # also 'bottom' and 'wake'; arrays ordered downstream
top['x'], top['ue'], top['dstar'], top['theta'], top['cf'], top['hk'], top['hstar']
top['n']                         # amplification factor on laminar stations (NaN when turbulent)
top['ctau']                      # max shear coefficient on turbulent stations (NaN when laminar)
nodes = xf.get_bl_nodes()        # the same per node in airfoil order + wake, Ue signed as in DUMP
```

Useful e.g. for the separation margin of a pressure recovery: the shape factor `hk` rises and `cf` drops
towards zero as a turbulent boundary layer approaches separation.

As in XFOIL, two airfoils are kept: the airfoil you set (the "buffer", returned by `xf.airfoil`) and the paneled
airfoil the analysis uses. `repanel()` builds the paneled airfoil from the buffer and never changes the buffer.

### Differences from upstream xfoil-python

These come from the reference version (`docs/xfoil-python.patch`), which xfoil_c reproduces:

- `a()` and `cl()` return 6 values: also the transition points `xtr_top`, `xtr_bot` (upstream: 4).
- `get_cp_distribution()` returns `x, y, cp` (upstream: `x, cp`).
- Double precision at the library boundary (upstream passes float32, so its results differ slightly).
- Setting `M` also sets XFOIL's `MINF1` (upstream does not).

### Differences from the reference

Deliberate, everything else is identical:

- `filter()` is removed: through the API it only changed inverse-design data that nothing uses.
- `aseq`/`cseq` give exactly the points of `range(start, end, step)`. The original could lose a point to rounding
  (`aseq(0, 0.7, 0.1)` gave 6 points with step 0.1167) and then stretched the step over the range.
- `get_cp_distribution()` returns NaN for cp unless the last analysed point converged and nothing changed the
  airfoil, paneling, flap or Mach number since (the original returned the stale or unconverged values).
- Every `XFoil()` instance is independent (see above).
- Inputs that made the Fortran code read or write outside its arrays now raise `XFoilError` or are handled
  (e.g. more than 2560 input points, a flap hinge behind the trailing edge).

## Build variants

The default build is exact. Faster variants are chosen at build time; the header of `CMakeLists.txt` documents
both switches in full.

```bash
pip install .                                                # exact (default)
XFOIL_ARCH=native XFOIL_FAST=ON pip install .                # fastest on this machine
python -c "import xfoil_c; print(xfoil_c.build_info())"      # which variant is installed
```

- `XFOIL_ARCH=native` (or any `-march` value such as `x86-64-v3`) uses your CPU's instructions including FMA.
  Results differ from the exact build in the last digits; a `native` build may not run on other CPUs.
- `XFOIL_FAST=ON` enables source-level fast paths (for example a vectorised boundary-layer solver).

## Tests

```bash
python setup.py build_ext --inplace        # build the library into xfoil_c/
python -m pytest tests -q
```

On its own the suite checks the library against recorded reference results (`tests/golden/`) and runs the
API, error-handling and thread-safety tests. Most tests compare every ported routine bit for bit with the
Fortran library; they need an [xfoil-python](https://github.com/daniel-de-vries/xfoil-python) checkout with
`docs/xfoil-python.patch` applied (its header says how), next to this repository (or
`XFOIL_PYTHON_DIR=/path/to/xfoil-python`), with the Fortran library built, and are skipped otherwise. The patch
also builds the Fortran library without auto-vectorisation (`-fno-tree-vectorize`): otherwise gfortran uses
glibc's vector math, which differs in the last bits.

The results also depend on the platform's math library (`exp`, `pow`, `atan2`, ... differ in the last bits
between glibc versions and the Windows C runtime). The golden file records the platform it was made on, and
the bit-for-bit comparison with it (`tests/test_regression.py`) runs only there. Everywhere else, and in a
non-exact build variant, `tests/test_reference.py` checks the results against it within tolerance and
requires no more unconverged points. It leaves out one sequence whose convergence depends on the last bits
(a `cseq` started from a stalled boundary layer, see `tests/golden/regression_cases.py`, case 1).

## Project layout

- `src/`: the C sources. `xfoil_c.h` declares the exported API; the solver state of one instance is the struct
  `xfoil_t` (`xstate.h`, generated from the Fortran modules by `tools/gen_globals.py`).
- `xfoil_c/`: the Python package. `xfoil.py` is the original xfoil-python wrapper with a marked extension
  section; `model.py` is unchanged; `_lib.py` loads the library with ctypes.
- `tests/`: pytest suite (see above).
- `tools/`: code generator for the state struct, Fortran call-graph tool, benchmark against xfoil-python.
- `docs/PORTING.md`: how the port was done, the design decisions and the list of XFOIL quirks that are
  reproduced on purpose; `docs/API_SNAPSHOT.md`: the xfoil-python API this port reproduces. The C API is
  documented in `src/xfoil_c.h`.

## Credits and license

XFOIL is by Mark Drela and Harold Youngren (MIT), https://web.mit.edu/drela/Public/web/xfoil/.
xfoil-python, the Fortran 90 module and Python wrapper this port is based on, is by Daniël de Vries,
https://github.com/daniel-de-vries/xfoil-python. The flap routine follows XFOIL 6.99's GDES FLAP.

The port was made with the help of Claude (Anthropic).

Licensed under the GNU General Public License v3 or later, as XFOIL and xfoil-python (see `LICENSE`).
