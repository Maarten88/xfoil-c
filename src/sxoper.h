/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of s_xoper.f90.
 */
#ifndef XFOIL_C_SXOPER_H
#define XFOIL_C_SXOPER_H

#include "xfoil_c.h"
#include "xstate.h"

/* Flap hinge moment hmom and forces hfx, hfy about (xof, yof) from the surface cp
 * (cpv if lvisc, else cpi). Requires lflap: the Fortran branch for lflap = .false.
 * asks for the hinge location on the console (getxyf) and is not ported; every
 * caller checks lflap first. */
XFOIL_API void mhinge(xfoil_t *xf);

#endif /* XFOIL_C_SXOPER_H */
