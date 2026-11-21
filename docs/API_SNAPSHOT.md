# API Snapshot: xfoil-python (contract for the C port)

Captured 2026-11-20 (Phase 0, Step 4) by introspecting `xfoil.xfoil` / `xfoil.model` and running each method
(gfortran 14.2, double-precision boundary). `xfoil_c` must reproduce this surface exactly.

## 1. Python package surface

`import xfoil` exposes `XFoil`, and the submodules `xfoil.xfoil` and `xfoil.model`.
`xfoil.xfoil` also has the module-level names `fptr = POINTER(c_double)` and `bptr = POINTER(c_bool)`.

### class `xfoil.model.Airfoil(x, y)`
| Member | Kind | Type / shape |
|---|---|---|
| `coords` | attribute | `ndarray (n, 2) float64` |
| `n_coords` | property (get only) | `int` (= `coords.shape[0]`) |
| `x`, `y` | property (get/set) | `ndarray (n,) float64` (column views of `coords`) |

### class `xfoil.xfoil.XFoil()`
Constructor: loads the library and calls `init()`. Destructor `__del__`: drops the handle (`FreeLibrary` on Windows).

**Properties** (all get/set)
| Property | Python type | Default after `init()` | Library calls (get / set) |
|---|---|---|---|
| `print` | `bool` | `True` | `get_print` / `set_print` |
| `airfoil` | `Airfoil` | none until set | `get_n_coords`+`get_airfoil` / `set_airfoil` |
| `Re` | `float` | (set to 1e6 in tests) | `get_reynolds` / `set_reynolds` |
| `M` | `float` | `0.0` | `get_mach` / `set_mach` |
| `xtr` | `tuple[float, float]` (top, bot) | - | `get_xtr` / `set_xtr` |
| `n_crit` | `float` | `9.0` | `get_n_crit` / `set_n_crit` |
| `max_iter` | `int` | `20` | `get_max_iter` / `set_max_iter` |

**Methods**
| Signature | Returns | Library call |
|---|---|---|
| `naca(specifier)` | `None` | `set_naca(int(specifier))` |
| `reset_bls()` | `None` | `reset_bls` |
| `repanel(n_nodes=160, cv_par=1, cte_ratio=0.15, ctr_ratio=0.2, xt_ref=(1, 1), xb_ref=(1, 1))` | `None` | `repanel` |
| `filter(factor=0.2)` | `None` | `filter` |
| `a(a)` | `tuple[6] float`: `cl, cd, cm, cp, xtr_top, xtr_bot` | `alfa` |
| `cl(cl)` | `tuple[6] float`: `a, cd, cm, cp, xtr_top, xtr_bot` | `cl` |
| `aseq(a_start, a_end, a_step)` | `tuple[5] ndarray (n,) float64`: `a, cl, cd, cm, cp` | `aseq` |
| `cseq(cl_start, cl_end, cl_step)` | `tuple[5] ndarray (n,) float64`: `a, cl, cd, cm, cp` | `cseq` |
| `get_cp_distribution()` | `tuple[3] ndarray (n_cp,) float64`: `x, y, cp` | `get_n_cp` + `get_cp` |

### Behavioural details the port must keep
- **Sequence length:** `n = abs(int((end - start) / step))`, computed in Python. The end point is excluded: `aseq(-20, 20, 0.5)` gives 80 points, and `cseq(-0.5, 0.5, 0.05)` gives 20.
- **Non-convergence:** `a`/`cl` return all six values as `nan`. In `aseq`/`cseq`, every array (including the swept `a`) is set to `nan` at points where `conv` is false. Example: `aseq(-20, 20, 0.5)` on NACA0012 has 2 `nan` points, including `a[0]`.
- **`naca`:** the specifier goes through `int()`, so leading zeros are lost (`'0012'` -> 12). Example: `naca('2412')` produces 425 coordinates.
- **Shared state:** `cdll.LoadLibrary` returns the same process-wide handle, so all `XFoil` instances share one global Fortran state.
- **Docstring mismatch (not the contract):** the docstring of `a()` lists 4 return values, but it actually returns 6.

## 2. C ABI (`libxfoil_fortran.so`)

Fortran `bind(c)` functions: every argument is passed **by pointer**, `real` = `double`, `logical(c_bool)` = `bool` (1 byte).
The exported symbols match `nm -D` exactly (27 in total).

```c
/* m_xfoil.f90 */
void   init(void);
void   xfoil(void);                       /* interactive program entry, not used by Python */

/* api.f90 */
void   set_print(const bool *setting);
bool   get_print(void);
void   set_airfoil(const double *x, const double *y, const int *n);
int    get_n_coords(void);
void   get_airfoil(double *x, double *y, const int *n);
void   set_naca(const int *spec);
void   set_reynolds(const double *Re);
double get_reynolds(void);
void   set_mach(const double *M);
double get_mach(void);
void   set_xtr(const double *xtr_top, const double *xtr_bot);
void   get_xtr(double *xtr_top, double *xtr_bot);
void   set_n_crit(const double *n_crit);
double get_n_crit(void);
void   set_max_iter(const int *max_iter);
int    get_max_iter(void);
void   reset_bls(void);
void   repanel(const int *n_panel, const double *cv_par, const double *cte_ratio, const double *ctr_ratio,
               const double *xs_ref1, const double *xs_ref2, const double *xp_ref1, const double *xp_ref2);
void   filter(const double *factor);
void   alfa(const double *a, double *cl, double *cd, double *cm, double *cp,
            double *xtr_top, double *xtr_bot, bool *conv);
void   cl  (const double *cl, double *a, double *cd, double *cm, double *cp,
            double *xtr_top, double *xtr_bot, bool *conv);
void   aseq(const double *a_start, const double *a_end, const int *n,
            double *a, double *cl, double *cd, double *cm, double *cp, bool *conv);   /* arrays length n */
void   cseq(const double *cl_start, const double *cl_end, const int *n,
            double *a, double *cl, double *cd, double *cm, double *cp, bool *conv);   /* arrays length n */
int    get_n_cp(void);
void   get_cp(double *x, double *y, double *cp, const int *n);
```

## 3. Reference values (NACA0012 from `xfoil/test.py`, Re=1e6, M=0, max_iter=40)
| Call | Result |
|---|---|
| `a(10)` | cl=1.0809, cd=0.0150, cm=0.0053, cp=-5.5077 |
| `cl(1)` | a=9.0617, cd=0.0135, cm=0.0013, cp=-4.7361 |
| `aseq(-20,20,0.5)`, `cseq(-0.5,0.5,0.05)` | arrays in `xfoil/test.py` (4 dp) |
