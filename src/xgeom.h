/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Airfoil geometry routines, ported from m_xgeom.f90 (the 7 routines that are
 *  reachable from the Python API; inter, interx, bendump, bendump2, ijsect and
 *  half are not ported).
 *  Conventions as in xspline.h: 0-based arrays and index values, double precision,
 *  scalar inputs by value, in/out scalars by pointer, Fortran argument order.
 *  Exported (XFOIL_API) so the tests can compare them with the Fortran originals.
 */
#ifndef XFOIL_C_XGEOM_H
#define XFOIL_C_XGEOM_H

#include "xfoil_c.h"

/* Leading-edge spline parameter *sle: where the surface tangent is normal to
 * the chord line to the TE midpoint. Splines x(s), y(s) given by x, xp, y, yp, s. */
XFOIL_API void lefind(double *sle, const double *x, const double *xp, const double *y,
                      const double *yp, const double *s, int n);

/* Arc length *sopp of the point opposite si on the other side of the chord line. */
XFOIL_API void sopps(double *sopp, double si, const double *x, const double *xp, const double *y,
                     const double *yp, const double *s, int n, double sle);

/* Splines the shape and scales it to unit chord with the LE at the origin
 * (x, y, s rescaled in place; xp, yp set). */
XFOIL_API void norm(double *x, double *xp, double *y, double *yp, double *s, int n);

/* Geometric parameters of the splined shape; t (n) is work space set to 1. */
XFOIL_API void geopar(const double *x, const double *xp, const double *y, const double *yp,
                      const double *s, int n, double *t, double *sle, double *chord, double *area,
                      double *radle, double *angte, double *ei11a, double *ei22a, double *apx1a,
                      double *apx2a, double *ei11t, double *ei22t, double *apx1t, double *apx2t,
                      double *thick, double *cambr);

/* Area/skin properties of the closed polygon x, y: itype 1 integrates over the
 * area, itype 2 over the skin with thickness t. */
XFOIL_API void aecalc(int n, const double *x, const double *y, const double *t, int itype,
                      double *area, double *xcen, double *ycen, double *ei11, double *ei22,
                      double *apx1, double *apx2);

/* Max thickness and camber evaluated at the discrete points. */
XFOIL_API void tccalc(const double *x, const double *xp, const double *y, const double *yp,
                      const double *s, int n, double *thick, double *xthick, double *cambr,
                      double *xcambr);

/* Max panel corner angle *amax (deg) at point *imax (0-based). iprint: 0 none,
 * 1 print the max, 2 also print every angle (behind show_output; printed
 * indices are 1-based, as in the Fortran). */
XFOIL_API void cang(const double *x, const double *y, int n, int iprint, int *imax, double *amax);

#endif /* XFOIL_C_XGEOM_H */
