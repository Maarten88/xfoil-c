/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Cubic spline routines, ported from m_spline.f90 (the API-reachable ones).
 *
 *  Conventions (whole port): arrays and index values are 0-based; reals are
 *  double; scalar inputs are passed by value and scalar in/out arguments by
 *  pointer. Argument order follows the Fortran routine.
 *  Exported (XFOIL_API) so the tests can compare them with the Fortran originals.
 */
#ifndef XFOIL_C_XSPLINE_H
#define XFOIL_C_XSPLINE_H

#include "xfoil_c.h"

/* As spline, with end conditions xs1/xs2: a specified dx/ds, 999.0 for zero 2nd
 * derivative, or -999.0 for zero 3rd derivative. n <= 1000. */
XFOIL_API void splind(const double *x, double *xs, const double *s, int n, double xs1, double xs2);

/* Solves the kk-long tridiagonal system (b sub-, a main, c super-diagonal).
 * d is replaced by the solution; a and c are destroyed. */
XFOIL_API void trisol(double *a, const double *b, double *c, double *d, int kk);

/* x(ss), dx/ds(ss) and d2x/ds2(ss) from spline arrays x, xs, s. */
XFOIL_API double seval(double ss, const double *x, const double *xs, const double *s, int n);
XFOIL_API double deval(double ss, const double *x, const double *xs, const double *s, int n);
XFOIL_API double d2val(double ss, const double *x, const double *xs, const double *s, int n);

/* Curvature of the splined 2-D curve (x(s), y(s)) at ss, and its derivative d(curv)/ds. */
XFOIL_API double curv(double ss, const double *x, const double *xs, const double *y,
                      const double *ys, const double *s, int n);

/* Inverse spline: Newton-iterates *si (initial guess in, result out) so that x(*si) = xi.
 * On failure prints a message (show_output) and restores the initial guess. */
XFOIL_API void sinvrt(double *si, double xi, const double *x, const double *xs, const double *s, int n);

/* Arc-length array s for the points (x, y). */
XFOIL_API void scalc(const double *x, const double *y, double *s, int n);

/* As spline / splind, but segments joined at identical successive s values
 * (derivative discontinuities allowed). Stops if the first or last point is duplicated. */
XFOIL_API void segspl(const double *x, double *xs, const double *s, int n);

#endif /* XFOIL_C_XSPLINE_H */
