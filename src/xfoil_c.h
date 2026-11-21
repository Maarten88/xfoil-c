/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Public C API of libxfoil_c. Exported functions use C linkage and, like the
 *  Fortran bind(c) API they replace, take all arguments by pointer.
 *  All reals are double precision.
 */
#ifndef XFOIL_C_H
#define XFOIL_C_H

#include <stdbool.h>

/* One XFoil instance: the whole solver state (src/xstate.h). */
typedef struct xfoil xfoil_t;

#if defined(_WIN32)
#  define XFOIL_API __declspec(dllexport)
#else
#  define XFOIL_API __attribute__((visibility("default")))
#endif

/* Library version string (matches xfoil_c.__version__). */
XFOIL_API const char *xfoil_version(void);
/* "exact" (bit-identical to the Fortran reference) or the build variant, e.g. "arch=native fast=ON" */
XFOIL_API const char *xfoil_build_info(void);

/* ---- api.f90 (contract: API_SNAPSHOT.md) ---- */
/* Contexts: every API function below takes the instance it works on as first argument. */
XFOIL_API xfoil_t *xfoil_new(void);           /* initial state, as Fortran before init(); NULL if out of memory */
XFOIL_API void xfoil_free(xfoil_t *xf);
XFOIL_API void init(xfoil_t *xf);             /* defaults (Fortran init) */

XFOIL_API void set_print(xfoil_t *xf, const bool *setting);
XFOIL_API bool get_print(xfoil_t *xf);
XFOIL_API void set_airfoil(xfoil_t *xf, const double *x_in, const double *y_in, const int *n_in);
XFOIL_API int get_n_coords(xfoil_t *xf);
XFOIL_API void get_airfoil(xfoil_t *xf, double *x_out, double *y_out, const int *n_out);
XFOIL_API void set_naca(xfoil_t *xf, const int *spec);
XFOIL_API void set_reynolds(xfoil_t *xf, const double *re);
XFOIL_API double get_reynolds(xfoil_t *xf);
XFOIL_API void set_mach(xfoil_t *xf, const double *m);
XFOIL_API double get_mach(xfoil_t *xf);
XFOIL_API void set_xtr(xfoil_t *xf, const double *xtr_top, const double *xtr_bot);
XFOIL_API void get_xtr(xfoil_t *xf, double *xtr_top, double *xtr_bot);
XFOIL_API void set_n_crit(xfoil_t *xf, const double *n_crit);
XFOIL_API double get_n_crit(xfoil_t *xf);
XFOIL_API void set_max_iter(xfoil_t *xf, const int *max_iter);
XFOIL_API int get_max_iter(xfoil_t *xf);
XFOIL_API void reset_bls(xfoil_t *xf);
XFOIL_API void repanel(xfoil_t *xf, const int *n_panel, const double *cv_par, const double *cte_ratio, const double *ctr_ratio,
                       const double *xs_ref1, const double *xs_ref2, const double *xp_ref1, const double *xp_ref2);
XFOIL_API void alfa(xfoil_t *xf, const double *a_input, double *cl_out, double *cd_out, double *cm_out, double *cp_out,
                    double *xtrtop_out, double *xtrbot_out, bool *conv);
XFOIL_API void cl(xfoil_t *xf, const double *cl_input, double *a_out, double *cd_out, double *cm_out, double *cp_out,
                  double *xtrtop_out, double *xtrbot_out, bool *conv);
XFOIL_API void aseq(xfoil_t *xf, const double *a_start, const double *a_end, const int *n_step, double *a_arr, double *cl_arr,
                    double *cd_arr, double *cm_arr, double *cp_arr, bool *conv_arr);
XFOIL_API void cseq(xfoil_t *xf, const double *cl_start, const double *cl_end, const int *n_step, double *a_arr,
                    double *cl_arr, double *cd_arr, double *cm_arr, double *cp_arr, bool *conv_arr);
XFOIL_API int get_n_cp(xfoil_t *xf);
XFOIL_API void get_cp(xfoil_t *xf, double *x_out, double *y_out, double *cp_out, const int *n_points);

/* Hard errors (Fortran `stop`): the public call returns early and the message is kept until
 * read. Returns 1 and copies the message (truncated to len) if an error occurred in the last
 * call, else 0; reading clears it. The Python package raises xfoil_c.XFoilError. */
XFOIL_API int xfoil_error(char *buf, int len);

/* ---- extensions (not in api.f90; Phase 10b) ---- */

/* Geometry of the current (paneled) airfoil, as XFoil's GEOPAR/TCCALC compute it; does not
 * change the state. thick, camber: max thickness and max camber, taken at the discrete nodes
 * (TCCALC), with x_thick, x_camber the x coordinate of the matching opposite-side point;
 * le_radius: LE radius (0 if the LE curvature is tiny); te_angle: TE angle in radians;
 * area: enclosed area. Lengths in airfoil coordinates (fractions of chord for a unit-chord
 * airfoil). All NaN if no airfoil is set (n < 2). */
XFOIL_API void get_geometry(xfoil_t *xf, double *thick, double *x_thick, double *camber, double *x_camber,
                            double *le_radius, double *te_angle, double *area);

/* Flap hinge (XFOIL GDES FLAP / OPER FMOM without prompts). set_hinge: hinge at x_hinge and
 * y_hinge (absolute, or y/t of the local buffer-airfoil thickness, 0 bottom .. 1 top, if
 * *y_relative) on the buffer and the current airfoil; hinge moments then follow every analysis.
 * No-op without airfoil; with a deflected flap, the deflection is re-applied about the new hinge.
 * flap: set the flap deflection to *angle degrees (+ = TE down; absolute: the undeflected airfoil,
 * kept at the first flap, is deflected once; 0 = undeflected), then repanel (pangen) with the
 * current parameters. No-op without hinge.
 * No-op without hinge. get_hinge: hinge moment and forces of the last analysed point (NaN without
 * hinge or before the next analysis). get_hinge_point: hinge x, y (NaN if none); returns 1 if a
 * hinge is set. clear_hinge: remove the hinge (and forget the undeflected airfoil). */
XFOIL_API void set_hinge(xfoil_t *xf, const double *x_hinge, const double *y_hinge, const bool *y_relative);
XFOIL_API void flap(xfoil_t *xf, const double *angle);
XFOIL_API void get_hinge(xfoil_t *xf, double *hinge_moment, double *hinge_fx, double *hinge_fy);
XFOIL_API int get_hinge_point(xfoil_t *xf, double *x_hinge, double *y_hinge);
XFOIL_API void clear_hinge(xfoil_t *xf);
/* current flap deflection in degrees (flap() sets it; 0 after a new airfoil) */
XFOIL_API double get_flap_angle(xfoil_t *xf);

#endif /* XFOIL_C_H */
