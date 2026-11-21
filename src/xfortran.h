/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Fortran intrinsic semantics needed for bit-identical results. Verified
 *  against gfortran 14 (-O2 -fdefault-real-8) with a probe program:
 *    MAX(a,b) == (a > b) ? a : b   (ties and +-0 give b; NaN only from b)
 *    MIN(a,b) == (a < b) ? a : b
 *    SIGN(a,b) == copysign(a, b)
 *  (Messages are plain printf since clean-up C2; no Fortran format emulation.)
 */
#ifndef XFOIL_C_XFORTRAN_H
#define XFOIL_C_XFORTRAN_H

#include <math.h>

static inline double fmax_f(double a, double b) { return (a > b) ? a : b; }
static inline double fmin_f(double a, double b) { return (a < b) ? a : b; }

#endif /* XFOIL_C_XFORTRAN_H */
