> This is the plan and decision record kept while porting xfoil-python's Fortran core to C. It was written for
> the work repository, where `xfoil-c/` sat next to the original `xfoil-python/`, so paths in it are relative to
> that layout; the per-step migration log it mentions is not part of this repository. The sections on the
> clean-up phase ("Phase 12") list every XFOIL quirk the port reproduces on purpose.

# XFOIL-Python → C-Backed Python Module: Migration Plan

## Overview

**Goal:** Port `xfoil-python/` (Fortran-backed) to `xfoil-c/` (C-backed) as a parallel Python module folder, sitting alongside the original.

**Source structure (Fortran):**
- `xfoil-python/src/` — ~25 Fortran modules (`.f90`): XFOIL internals (panel method, boundary layer, geometry, design, solver, I/O, etc.)
- `xfoil-python/xfoil/` — Python wrapper (`xfoil.py`, `model.py`, `__init__.py`, `test.py`)
- Build: `setup.py` + `CMakeLists.txt` + `pyproject.toml`

**Target structure (C):**
- `xfoil-c/src/` — C translation of each Fortran module (`.c` / `.h`)
- `xfoil-c/xfoil_c/` — Python package mirroring `xfoil/` API surface
- `xfoil-c/` — `setup.py` / `pyproject.toml` / `CMakeLists.txt` for C extension build

**Key constraints:**
- No Fortran compiler installed; Phase 0 must also verify/repair the original Fortran build as a baseline.
- Each phase must fit in one short working session (~30–90 min) and end with an update to `MIGRATION_LOG.md`.
- Do NOT read any Fortran source files until the relevant phase.
- Keep the public Python API identical so that downstream code can swap `import xfoil` → `import xfoil_c`.

## Scope update (2026-11-20, user decision)

Only routines reachable from the Python API are ported (121 of 208 routines, ~11.8k of 22k lines), as listed by `tools/fortran_callgraph.py`. This removes **m_aread** (Phase 2), **m_xpol** (Phase 5) and **m_iopol** (Phase 7) entirely, and most of m_xmdes, m_xqdes (Phase 8), m_xoper and m_xfoil (Phases 7 and 9). Interactive prompts (`aski`/`askr`) are replaced by non-interactive behaviour; see CLAUDE.md.

---

## Phase 0 — Build & Test Baseline

**Session:** 1  
**Goal:** Get the existing Fortran project to compile and run, establishing a working reference.

### Steps
1. Install a Fortran toolchain (e.g. `gfortran`, `cmake`, and `ninja-build` via `apt`).
2. Attempt the existing build:
      cd xfoil-python
      pip install .
   Debug CMake / build errors. Save the working build config.
3. Test if `xfoil/test.py` loads the xfoil module succesfully, fix the python test code or module if necessary. 
   Then run `xfoil/test.py` against the reference data files in `runs/` and record which tests pass/fail.
4. Capture a snapshot of the public Python API (function names, signatures, return shapes) by introspecting `xfoil.xfoil`.
5. Create `xfoil-c/` directory skeleton:
   ```
   xfoil-c/
   ├── pyproject.toml
   ├── setup.py
   ├── CMakeLists.txt
   ├── README.md
   ├── xfoil_c/
   │   ├── __init__.py
   │   └── model.py          (copy from xfoil/, will adapt later)
   └── src/
   ```
6. Create `MIGRATION_LOG.md` with initial entry.

**Exit criteria:** Fortran build succeeds, tests run, baseline API captured, `xfoil-c/` scaffold exists.

---

## Phase 1 — C Extension Skeleton (pybind11 / cffi)

**Session:** 2  
**Goal:** Wire up a trivial C extension that imports into `xfoil_c` and proves the build toolchain.

### Steps
1. Add `pybind11` (or `cffi`) to `xfoil-c/pyproject.toml` build-requires.
2. Write a minimal `src/xfoil_stub.c` + `src/xfoil_stub.h` exporting one function `xfoil_version()` that returns a string.
3. Configure `setup.py` / `CMakeLists.txt` to compile it as a Python extension module `_xfoil_c`.
4. In `xfoil_c/__init__.py`, import `_xfoil_c` and expose the stub.
5. Write `xfoil-c/tests/test_skeleton.py` asserting the stub import works.
6. Update `MIGRATION_LOG.md`.

**Exit criteria:** `pip install ./xfoil-c` works; `import xfoil_c` succeeds; stub test passes.

---

## Phase 2 — Utility & Low-Level Modules (no airfoil logic)

**Session:** 3  
**Modules to port:**
- `m_sort.f90` → `xsort.c` / `xsort.h`
- `m_spline.f90` → `xspline.c` / `xspline.h`
- `m_xutils.f90` → `xutils.c` / `xutils.h`
- `m_userio.f90` → `xuserio.c` / `xuserio.h`
- `m_aread.f90` → `xaread.c` / `xaread.h`

These are self-contained helpers (sorting, interpolation, file I/O, utility routines) with minimal coupling.

### Steps
1. For each module: read the `.f90`, translate to C, compile as part of the extension.
2. Bind each C function into the Python extension layer.
3. Unit-test each translated module against its Fortran counterpart output (write small C test harnesses or compare via Python).
4. Update `MIGRATION_LOG.md`.

**Exit criteria:** 5 utility modules compiled, bound, and validated against Fortran reference.

---

## Phase 3 — Geometry & NACA Generation

**Session:** 4  
**Modules to port:**
- `m_naca.f90` → `xnaca.c` / `xnaca.h`
- `m_xgeom.f90` → `xgeom.c` / `xgeom.h`
- `i_circle.f90` → `xcircle.c` / `xcircle.h`
- `i_pindex.f90` → `xpindex.c` / `xpindex.h`

### Steps
1. Translate geometry routines (airfoil shape generation, circle mapping, panel indexing).
2. Bind into Python extension.
3. Validate: generate NACA 4-digit airfoils in C and diff against Fortran output.
4. Update `MIGRATION_LOG.md`.

**Exit criteria:** NACA/geometry functions ported and validated.

---

## Phase 4 — Panel Method

**Session:** 5  
**Modules to port:**
- `m_xpanel.f90` → `xpanel.c` / `xpanel.h`
- Related: portions of `m_xgeom.f90` if panel-specific logic lives there.

### Steps
1. Port panel discretisation, influence coefficients, singularity distribution.
2. Bind and test: generate panel distributions for a known airfoil, compare to Fortran.
3. Update `MIGRATION_LOG.md`.

**Exit criteria:** Panel method C code compiles and produces matching panel data.

---

## Phase 5 — Solver Core (Linear System)

**Session:** 6  
**Modules to port:**
- `m_xsolve.f90` → `xsolve.c` / `xsolve.h`
- `m_xpol.f90` → `xpol.c` / `xpol.h`

### Steps
1. Port the linear solver / polynomial evaluation routines.
2. These are mathematically dense; verify numerics with small matrix test cases.
3. Update `MIGRATION_LOG.md`.

**Exit criteria:** Solver routines compiled and numerically validated.

---

## Phase 6 — Boundary Layer System

**Session:** 7 (may split into 6a / 6b if too large)  
**Modules to port:**
- `m_xbl.f90` → `xbl.c` / `xbl.h`
- `m_xblsys.f90` → `xblsys.c` / `xblsys.h`
- `i_blpar.f90` → `xblpar.c` / `xblpar.h`
- `i_xbl.f90` → `xibln.c` / `xibln.h`
- `s_xbl.f90` → `sxbl.c` / `sxbl.h`

### Steps
1. Port the Thwaites / epr boundary-layer formulation and system assembly.
2. Validate against Fortran BL output on a known airfoil at known α.
3. Update `MIGRATION_LOG.md`.

**Exit criteria:** BL modules compile; BL coefficients match Fortran within tolerance.

---

## Phase 7 — Operating Modes & Polars

**Session:** 8  
**Modules to port:**
- `m_xoper.f90` → `xoper.c` / `xoper.h`
- `m_iopol.f90` → `xiopol.c` / `xiopol.h`
- `s_xoper.f90` → `sxoper.c` / `sxoper.h`

### Steps
1. Port operating-mode logic (direct/inverse modes, α/force specification).
2. Port polar I/O and polar file formatting.
3. Validate: run a single α sweep and compare lift/drag/moment to Fortran.
4. Update `MIGRATION_LOG.md`.

**Exit criteria:** Operating modes and polar I/O ported and validated.

---

## Phase 8 — Design Modules

**Session:** 9  
**Modules to port:**
- `m_xgdes.f90` → `xgdes.c` / `xgdes.h`
- `m_xmdes.f90` → `xmdes.c` / `xmdes.h`
- `m_xqdes.f90` → `xqdes.c` / `xqdes.h`

### Steps
1. Port geometry-design, mode-design, and quantity-design routines.
2. Validate: perform a single design iteration (e.g. inverse Cp target → geometry update) and compare to Fortran.
3. Update `MIGRATION_LOG.md`.

**Exit criteria:** All design modules compiled and validated.

---

## Phase 9 — Main Orchestration & Entry Point

**Session:** 10  
**Modules to port:**
- `m_xfoil.f90` → `xfoil_core.c` / `xfoil_core.h`
- `p_xfoil.f90` → `xfoil_entry.c` / `xfoil_entry.h`
- `s_xfoil.f90` → `sxfoil.c` / `sxfoil.h`
- `api.f90` → `xfoil_api.c` / `xfoil_api.h`
- `i_xfoil.f90` → `ixfoil.c` / `ixfoil.h`

### Steps
1. Port the top-level XFOIL state machine, API dispatch, and entry points.
2. Wire all previously-ported modules into the C orchestrator.
3. Bind the full public API into the Python extension.
4. Run the full `xfoil/test.py` equivalent (`xfoil-c/tests/test_full.py`) and compare results.
5. Update `MIGRATION_LOG.md`.

**Exit criteria:** Full C-backed XFOIL runs end-to-end; results match Fortran within numerical tolerance.

---

## Phase 10 — Python API Alignment & Integration Tests

**Session:** 11  
**Goal:** Make `xfoil_c` a drop-in replacement for `xfoil`.

### Steps
1. Align `xfoil_c/__init__.py`, `xfoil_c/model.py`, and the extension module's public surface with the original `xfoil/` package.
2. Port `xfoil/test.py` → `xfoil-c/tests/test_xfoil_c.py`, running the same assertions against the same `runs/` data files.
3. Add a cross-validation test that runs both `xfoil` and `xfoil_c` on the same inputs and diffs results.
4. Update `MIGRATION_LOG.md`.

**Exit criteria:** `import xfoil_c` behaves identically to `import xfoil` for all public functions; all tests pass.

---

## Phase 10b - Extensions (user decision 2026-11-20; after the drop-in API is complete and verified)

New C API functions + `XFoil` members in `xfoil_c` only (the drop-in API stays identical; these are additions):
1. **Geometry properties (read-only):** max thickness and its x, max camber and its x (`tccalc`), LE radius,
   TE angle, area (`geopar`) of the current airfoil. Tested against the Fortran `tccalc`/`geopar`.
2. **Flap + hinge moment (user decisions 2026-11-20; steps 2+3 merged):**
   - Revised (user, 2026-11-20): ONE hinge. `xf.set_hinge(x_hinge, y_hinge, y_relative=False)` sets it on the buffer
     (XBF/YBF, LBFLAP) and on the current airfoil (XOF/YOF, LFLAP) - the non-interactive XFOIL FMOM/GETXYF path - so
     hinge moments work at 0 deg without repaneling; y_relative: y/t (0 bottom, 1 top) as GETXYF's "999" option.
   - `xf.flap(angle)`: literal XFOIL FLAP about that hinge (+ = TE down, degrees, cumulative, 0 = no-op, error without
     hinge), then automatic repanel with the current paneling parameters (pangen), so analysis and xf.airfoil agree.
   - Results: read-only `hinge_moment`, `hinge_fx`, `hinge_fy` of the last analysed point (mhinge inside
     specal/speccl/viscal when lflap); NaN while no hinge is set and until the next analysis after set_hinge/flap.
   - A new airfoil (`xf.airfoil = ...`, `xf.naca(...)`) clears the hinge (Python extension; Fortran would keep lflap).
   - Port from the upstream XFOIL 6.99 source (xgdes.f FLAP, INSIDE, GETXYF without prompts; xgeom.f SSS,
     SCHECK); plotting dropped. Reference for tests: those upstream routines extracted and compiled with the
     reference flags (-O2 -fdefault-real-8 -fno-tree-vectorize) plus plot stubs, compared bit for bit on the
     deflected buffer.
3. **Boundary layer (2026-11-21, after the move to this repository):** `get_n_bl`/`get_bl` (C) return per node
   (airfoil order, then wake) the quantities of m_xoper bldump (OPER DUMP): s, x, y, Ue/Vinf, dstar, theta, Cf,
   Hk, H*, plus N on laminar and Ctau on turbulent stations (both live in CTAU). NaN unless the last point is
   a converged viscous one. Python: `XFoil.get_bl_nodes()` (per node) and `XFoil.get_bl()` (split into top,
   bottom, wake, ordered downstream, |Ue|; None without a converged viscous point). Tested against the file
   the Fortran bldump writes after the same API calls (within its print precision; `tests/test_bl.py`).
Not wanted: surface smoothing (XFOIL only offers it through inverse design), inverse design itself.

## Phase 11 — Cleanup, Performance, and Polish

**Session:** 12  
### Steps
1. Remove any Fortran remnants from the C build.
2. Benchmark: compare wall-clock time of `xfoil` vs `xfoil_c` on representative workloads.
3. Profile hotspots; consider OpenMP threading for panel/solver loops if beneficial.
4. Finalise `pyproject.toml`, `setup.py`, `CMakeLists.txt` for clean `pip install ./xfoil-c`.
5. Update `MIGRATION_LOG.md` with final summary.

**Exit criteria:** Clean C-only build, performance report, no Fortran dependency.

---

## File Convention (C side)

| Fortran | C | Notes |
|---------|---|-------|
| `m_xbl.f90` | `src/xbl.c` + `src/xbl.h` | Module → header/impl pair |
| `i_xfoil.f90` | `src/ixfoil.c` + `src/ixfoil.h` | Interface file |
| `s_xbl.f90` | `src/sxbl.c` + `src/sxbl.h` | Subroutine-level helper |
| `p_xfoil.f90` | `src/xfoil_entry.c` + `src/xfoil_entry.h` | Program entry |

Naming: drop `m_`/`i_`/`s_`/`p_` prefixes; use lower-case snake-style where natural.

---

## Tracking: MIGRATION_LOG.md

Each phase update appends a section:

```markdown
## Phase N — <title>  (<date>)
- **Done:** bullet list of what was completed.
- **Modules touched:** list.
- **Build status:** pass / fail notes.
- **Open issues:** anything deferred.
- **Next phase preview:** what comes next.
```

---

## Dependency Graph (rough)

```
Phase 0 (baseline)
  └─ Phase 1 (C skeleton)
       └─ Phase 2 (utils: sort, spline, utils, userio, aread)
            ├─ Phase 3 (geometry: naca, geom, circle, pindex)
            │    └─ Phase 4 (panel: xpanel)
            │         └─ Phase 5 (solver: xsolve, xpol)
            │              └─ Phase 6 (boundary layer: xbl, xblsys, blpar, ibl, sxbl)
            │                   └─ Phase 7 (operating modes: xoper, iopol, sxoper)
            │                        ├─ Phase 8 (design: xgdes, xmdes, xqdes)
            │                        └─ Phase 9 (main orchestration)
            │                             └─ Phase 10 (API alignment + integration)
            │                                  └─ Phase 11 (cleanup + perf)
```

Each phase depends on the previous; phases 3 & 4 are sequential (geometry → panel). Phases 5–7 are strictly sequential (solver → BL → ops). Phase 8 can overlap with Phase 7 in a separate branch.

---

## Phase 12 - Clean-up and optimization (after a fully working, verified C port)

**User plan (2026-11-20):** once `xfoil_c` is a working drop-in equivalent, clean up everything that only exists because of the Fortran legacy, poor coding or carried-over bugs, then optimize (CUDA or other options). The bit-identical test suite against the Fortran library is the safety net for the clean-up: first make changes that keep results identical, then deliberately accept small numerical differences where a fix or optimization requires it.

### Clean-up steps (Phases 11 + 12 merged; one commit per step; user go-ahead 2026-11-20)

Rule: steps C1-C6 keep every public result bit-identical (checked by the C0 regression net and the live
comparison with the Fortran package); only steps marked (numerics) may change results, each with a documented
reason and tolerance. Internal whole-state tests (test_x*.py against the Fortran module routines) are kept as long
as the internals they check still exist; a step that removes or restructures those internals retires the affected
tests explicitly (listed in its log entry), since the API-level net covers the behaviour.

- **C0 Safety net + baseline:** API-level golden regression data (all public calls incl. extensions on several
  airfoils/conditions, exact bytes, NaN pattern) generated from the verified build and checked in; benchmark script
  (Fortran package vs xfoil_c) with baseline timings.
- **C1 Dead code:** Phase-1 stubs (xfoil_echo_double), unused ported routines and their tests, test-only exports moved
  out of the public symbol set, unused includes.
- **C2 Output/messages:** drop Fortran formatting emulation (xfortran.h fort_put_*), plain messages behind show_output;
  fatal `xfoil_stop` -> error codes surfaced to Python (no process exit).
- **C3 Unused state (done 2026-11-20):** remove globals the library never reads (plot/menu/polar-file variables set by init, ISPARS/MSES
  domain parameters, ...), i_circle/filter internals reduced to the observable behaviour.
- **C4 Readability:** Fortran 1-based accessor macros and SPAG-style control flow -> plain 0-based C; remove copies
  only needed for EQUIVALENCE (update), pointer-alias emulation (blcom_t) where it can be simplified.
  Sub-steps: C4a BL Newton-system blocks (VS1/VS2/VSREZ/VSM/VSR/VSX/VTMP, VA/VB/VM/VDEL/VZ in xblsys.c,
  xbl.c, xsolve.c) -> 0-based (done 2026-11-20); C4b geometry accessors (X_/Y_/S_/CP_, XB/YB/SB/W1..W5/SNEW,
  XS/YS in xfoil.c pangen, sxoper.c, xgdes.c, xfoil_ext.c) -> 0-based (done 2026-11-20); C4c goto/SPAG control flow (xbl, xblsys,
  xpanel, xfoil_ext) (done 2026-11-20); C4d update's EQUIVALENCE copies, blcom_t simplification (done 2026-11-20).
- **C5 State container:** module globals -> one context struct (prerequisite for multiple XFoil instances, threads,
  GPU); fixed array sizes -> sizes chosen at run time with error returns instead of `stop`.
  Sub-steps (user decisions 2026-11-20: explicit `xfoil_t *xf` parameter; each Python XFoil() gets independent
  state): C5a struct + xf parameter, API on a default context (done 2026-11-20); C5b API takes a context, Python
  XFoil() owns one, error/guard state per thread (done 2026-11-20); C5c run-time array sizes (allocation per context, error returns
  instead of NVX/IQX stops); C5d scratch statics -> context (xscratch.h), thread-safe default context, threads test (done 2026-11-20;
  kept: ~30 one-function Fortran module variables stay context fields (already per instance; they are compared
  with Fortran), and the test-only exports stay while the per-routine Fortran comparison tests exist).
- **C6 (numerics) Legacy bugs:** out-of-bounds reads, uninitialised values, unbounded rlx, ... fixed; results may
  change in edge cases; each fix documented.
  Sub-steps: C6a memory safety (done 2026-11-20); C6b uninitialised values (none found), dead computations (done 2026-11-20); C6c numerical/API
  quirks (rlx bounds, get_cp after unconverged points, aseq/cseq point count; done 2026-11-20, user decisions).
- **C7 (numerics) Build/arithmetic:** drop bit-identity-only constructs (`-(a*b)`, powi, MIN/MAX order, -ffp-contract=off,
  Fortran reference built with -fno-tree-vectorize) where they cost speed.
  Measured 2026-11-20 (6 flag variants, golden deviation + timing): -O3 -fno-math-errno bit-identical and adopted;
  FMA/-march=x86-64-v3/native 10-15% faster but up to 1.2e-5 relative deviation (and -O3 + v3: 250 convergence
  changes on the golden cases), needs AVX2: not adopted (user decision pending). powi / -(a*b) / MIN-MAX order /
  -ffp-contract=off cost nothing measurable (-O3 with them is bit-identical and as fast): kept. Profile (gprof,
  viscous aseq): blsolv 61%, setbl 7%, baksub 6%, source_half_panels 5%, qdcalc 5% -> optimization phase.
- **Then optimization:** profile, algorithmic improvements, OpenMP for sequences, CUDA evaluation.
  User decisions (2026-11-20): small deviations from the reference are fine and individual convergence outcomes may
  differ, but non-convergence must not become more frequent; keep the bit-identical exact build as the default and
  baseline; optimizations are build options (aggressive strategies welcome). Steps: O1 build variants + test
  framework (done 2026-11-20); O2 blsolv (61 % of a viscous sequence; done 2026-11-21: XFOIL_FAST row-contiguous VM, blsolv 3.5x,
  native+fast analyses at 0.66-0.68 of Fortran); next profile: blsolv 32 %, setbl 14 %, qdcalc 8 %, gauss 6 %,
  source_half_panels 6 %, baksub 5 %, update 5 %.
  O3 (2026-11-21): line profile shows the vectorised blsolv row update near the hardware limit (38 %); other
  lines < 5 %; fusion/helper tries gave 1-2 % and were reverted. Next: measure the real optimizer workload.

### Legacy quirks currently reproduced on purpose (clean-up candidates)
- (removed in C1, unreachable from the API: m_sort, spline, splina, splnxy with its `sqrt(sint)` bug, segspld, curvs)
- (fixed in C6a) `cang`: read `X(0)`/`X(N+1)` for doubled points next to the ends; such a point now has no corner (only the warning used it).
- (removed in C6b) Unused computations kept for a literal translation: `ggcalc` `psiinf`, `psilin` `qtanm`, `bldif` `uq_t1..uq_re` (+ `uq_upw`, `uq_rta`, `hr_rta`, `hkc_rta`), `trchek2` `z_xf` (+ `ax_xf`, `tt_xf`, `dt_xf`, `ut_xf`), `pangen` dead stores `cv1`/`cvs1`.
- `naca()` passes the designation through `int()` (leading zeros lost); `aski`/`askr` console prompts replaced by non-interactive behaviour.
- Fixed-size work arrays and hard limits with `stop` (since C2 an XFoilError in Python instead of a process exit) (NMAX = 1000 in splind, NVX = 500 in ludcmp, IVX/IZX/IQX overflow checks) -> dynamic sizes or error returns instead of terminating the process.
- Global state everywhere (i_xfoil, i_xbl, i_blpar): ~31 MB of statics, prevents several XFoil instances / threads; candidates for a context struct.
- COM1/COM2 (`blcom_t` union of a double array and named fields, i.e. the Fortran pointer aliases): kept, since the field names are the readable form and the array view gives the Fortran layout for the tests; since C4d copies are struct assignments (com1 = com2, c2sav = com2; c1sav/c2sav are blcom_t too).
- (removed in C4d) `update` (m_xbl): Fortran places its scratch arrays unew/u_ac/qnew/q_ac on copies of VA/VB (EQUIVALENCE) and copies them back, leaving VA/VB partly overwritten; C uses plain static locals and leaves VA/VB alone. Tests skip VA/VB in globals_diff (conftest.FORTRAN_SCRATCH) and compare their active blocks after setbl/blsolv (va_vb_diff).
- `mhinge` (s_xoper): Fortran reads `S(0)` (no panel with `S(i-1) < bots`) and, fixed in C6a, `S(N+1)` (no node past `tops`) out of bounds; C clamps to the first / last panel. Its `lflap = .false.` branch (getxyf console prompt) is not ported; all callers check lflap first.
- `viscal`: separation/reattachment summary (hkmax, psep, patt) is computed only for a printed line; C computes it only when printing.
- (removed in C3) `init`'s ~100 plotting/menu/polar-file variables, ISPARS/MSES domain parameters (and set_airfoil's uninitialised `Itype`), write-only results (psio, sigte_a/gamte_a, tindex, xssitr, yoctr, uslp, xcpmni/xcpmnv, dwte, gambl, ...): 231 globals, list in tools/dropped_globals.txt.
- `pangen` (m_xfoil): (fixed in C6a) read `SB(0)` if the LE lies between buffer points 1 and 2; now the TE equation at node 1 stays fixed there. `ipfac = 3` then `ipfac = 5` (C: 5); dsmin/dsmax computed for a commented-out print (dropped in C); `stop` when corners push N above IQX-1.
- api `filter` (Fortran): never initialises the circle-plane data (eiwset/scinit/cncalc are unreachable from the API), so it only sets cn(0) = NaN and the lqspec/lqsppl flags, which nothing reads. Since C3 the C `filter` was a no-op (i_circle, cnfilt, piqsum removed); then removed from the C library and the Python API (user decision, 2026-11-20; the only intended API difference).
- api `set_airfoil`: (fixed in C6a) no check of n_in against the buffer size IBX (memory overwrite above 2560 points): now XFoilError, also for fewer than 3 points. `set_naca(<=0)` writes to stderr even with printing off.
- api `aseq`/`cseq`: (fixed in C6c, user decision) the Python count abs(int((end-start)/step)) lost a point to rounding and the library then stretched the step over (end-start); now `_range_count` (ceil with tolerance) and the matching end, so the points are exactly range(start, end, step). The convergence flag mixes the caller's conv_arr value with lvconv; `get_cp` ignores n_points and writes N values (harmless: Python passes N); it computes cp from GAM with the Karman-Tsien formula, which equals CPV after a viscous run and CPI after an inviscid one; (fixed in C6c, user decision) it returned the last analysed point even if unconverged, and stale values after an airfoil/Mach change: now NaN cp unless the last point converged and nothing changed since (`xf->scr.cp_valid`). `alfa`/`cl` run the specal/speccl before the wake/convergence checks, `aseq` after. The exported interactive program `xfoil` (m_xfoil) is not ported.
- `flap` (ported from upstream XFOIL 6.99 FLAP): (fixed in C6a) break indices IT1/IT2/IB1/IB2 stayed undefined if no node brackets a break, then SB(IB1-1) etc. were read outside the buffer: now XFoilError (e.g. hinge behind the TE); SCHECK read X(N+1) while shifting (loop bound fixed). Buffer arrays keep stale values beyond NB after the shifts (harmless). Hinge state: XFoil keeps a buffer hinge (XBF/YBF, LBFLAP) and a current hinge (XOF/YOF, LFLAP); the extension API keeps them equal. Since 2026-11-21 (user decision) the API flap() is absolute: the undeflected buffer airfoil is kept in xf->scr (flap_base_*) at the first flap and restored before each deflection; set_hinge re-applies the current angle about a new hinge; flap_angle reports it.
- (fixed in C6c, user decision) `update`: with Re(CL) (MATYP 2/3) and CL < 0, DCLMIN = max(-0.5, -0.9 CL) is positive and the DCLMIN limit gave a negative or huge RLX; the alpha/CL limits now only reduce RLX to a positive value (normal cases unchanged; unreachable through the API, which keeps MATYP = 1).
- Index conventions: 0-based C vs 1-based Fortran mapping in tests; `tools/index_globals.txt`. `setbl` keeps the real global `tindex` (fractional transition panel index) on the Fortran 1-based scale.
- Unary minus written `-(a*b)`, `__builtin_powi` for integer powers, mirrored MIN/MAX operand order (mrchue `min(cti, 0.30)`): only needed for bit-identity with gfortran.
- (done in C2) Fortran output formatting helpers removed; messages are plain printf.
- Build: Fortran reference compiled with `-fno-tree-vectorize`; C with `-ffp-contract=off` (bit-identity only).
