/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Miscellaneous utilities, ported from m_xutils.f90.
 *  Conventions as in xspline.h: 0-based arrays, double precision, scalar inputs
 *  by value, in/out scalars by pointer, Fortran argument order.
 *  Exported (XFOIL_API) so the tests can compare them with the Fortran originals.
 */
#ifndef XFOIL_C_XUTILS_H
#define XFOIL_C_XUTILS_H

#include "xfoil_c.h"

/* Geometrically stretched array s[0..nn-1]: s[0] = 0, s[1] - s[0] = ds1,
 * s[nn-1] = smax, each increment `ratio` times the previous one.
 * Stops if nn < 3; prints a message (show_output) if the ratio does not converge. */
XFOIL_API void setexp(double *s, double ds1, double smax, int nn);

/* atan2(y, x) moved by a multiple of 2*pi to lie within pi of thold, so that
 * an angle array does not jump by 2*pi where the path crosses the -x axis. */
XFOIL_API double atanc(double y, double x, double thold);

#endif /* XFOIL_C_XUTILS_H */
