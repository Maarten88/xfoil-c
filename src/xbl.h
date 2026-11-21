/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  BL marching and Newton update, ported from m_xbl.f90 (API-reachable routines;
 *  added incrementally). preptrs is not needed: the COM1/COM2 aliases are the
 *  blcom_t union (ixbl.h). Side arguments `is` are 0-based (0 top, 1 bottom).
 *  Exported (XFOIL_API) so the tests can compare them with the Fortran originals.
 */
#ifndef XFOIL_C_XBL_H
#define XFOIL_C_XBL_H

#include "xfoil_c.h"
#include "xstate.h"

/* Limits *dstr so that the kinematic shape parameter is at least hklim. */
XFOIL_API void dslim(double *dstr, double thet, double uedg, double msq, double hklim);

/* BL parameter defaults (i_blpar). */
XFOIL_API void blpini(xfoil_t *xf);

/* Forced-transition BL coordinate xiforc for side is (from xstrip). */
XFOIL_API void xifset(xfoil_t *xf, int is);

/* Marches the BLs and wake in direct mode using uedg (inverse mode with an
 * extrapolated Hk where separation would occur), checking for transition.
 * Returns 0 if a NaN made the transition check fail (ABORT_ON_NAN), else 1. */
XFOIL_API int mrchue(xfoil_t *xf);

/* Marches the BLs and wake in mixed mode using the current Ue and Hk (along a line
 * quasi-normal to the Ue-Hk characteristic, avoiding the Goldstein/Levy-Lees
 * singularity), checking for transition. Returns 0 on a NaN abort, else 1. */
XFOIL_API int mrchdu(xfoil_t *xf);

/* Sets up the BL Newton system (va, vb, vm, vz, vdel) for the current BL variables and
 * edge velocities (initializing the BL with mrchue if needed, then mrchdu). Returns 0
 * on a NaN abort, else the last march / transition-check result. */
XFOIL_API int setbl(xfoil_t *xf);

/* Adds the Newton deltas (vdel) to the BL variables with under-relaxation, updates CL
 * (if lalfa) or alpha, and sets the rms/max changes (rmsbl, rmxbl, vmxbl, imxbl, ismxbl). */
XFOIL_API void update(xfoil_t *xf);

#endif /* XFOIL_C_XBL_H */
