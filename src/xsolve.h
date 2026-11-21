/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Linear solvers, ported from m_xsolve.f90 (ludcmp, baksub so far).
 *  Matrices use the Fortran column-major layout with leading dimension nsiz:
 *  Fortran A(i,j) is a[(j-1)*nsiz + (i-1)], so the i_xfoil globals (e.g. aij)
 *  can be passed directly as &aij[0][0].
 *  Pivot indices (indx) are 0-based; ludcmp and baksub use them consistently.
 *  Exported (XFOIL_API) so the tests can compare them with the Fortran originals.
 */
#ifndef XFOIL_C_XSOLVE_H
#define XFOIL_C_XSOLVE_H

#include "xfoil_c.h"
#include "xstate.h"

/* LU-factors the n x n matrix a in place (Crout, partial pivoting; Numerical
 * Recipes). indx receives the 0-based pivot rows. Stops if n > 500. */
XFOIL_API void ludcmp(int nsiz, int n, double *a, int *indx);

/* Solves A x = b with the factors from ludcmp; b is replaced by x. */
XFOIL_API void baksub(int nsiz, int n, const double *a, const int *indx, double *b);

/* Solves the nn x nn system z (leading dimension nsiz) for nrhs right-hand sides r
 * (column-major, leading dimension nsiz), replaced by the solutions; z is destroyed.
 * Gaussian elimination with partial pivoting. */
XFOIL_API void gauss(int nsiz, int nn, double *z, double *r, int nrhs);

/* Custom block solver for the coupled viscous-inviscid Newton system in the i_xfoil
 * globals va, vb, vm, vz, vdel (nsys lines); the solution replaces vdel. */
XFOIL_API void blsolv(xfoil_t *xf);

#endif /* XFOIL_C_XSOLVE_H */
