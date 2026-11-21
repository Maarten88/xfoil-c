/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of the API-reachable part of m_xoper.f90.
 */
#ifndef XFOIL_C_XOPER_H
#define XFOIL_C_XOPER_H

#include "xfoil_c.h"
#include "xstate.h"

/* Minimum cp (cpmni, cpmnv, cpmn) and its x location over airfoil and wake. */
XFOIL_API void fcpmin(xfoil_t *xf);

/* Inviscid solution at the specified alpha (Newton loop on M(CL) for matyp > 1). */
XFOIL_API void specal(xfoil_t *xf);

/* Inviscid solution at the specified CL clspec (Newton loop on alpha). */
XFOIL_API void speccl(xfoil_t *xf);

/* Converges the viscous operating point with up to niter1 Newton iterations
 * (niter1 <= 0: no iterations; Fortran prompts for the count when it is 0).
 * Returns 0 if setbl aborted on a NaN (ABORT_ON_NAN), else 1. Convergence is lvconv. */
XFOIL_API int viscal(xfoil_t *xf, int niter1);

#endif /* XFOIL_C_XOPER_H */
