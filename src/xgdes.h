/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of the API-reachable part of m_xgdes.f90 (getxyf is not ported: console prompts).
 */
#ifndef XFOIL_C_XGDES_H
#define XFOIL_C_XGDES_H

#include "xfoil_c.h"
#include "xstate.h"

/* Current airfoil (x, y, n) := buffer airfoil (xb, yb, nb) without doubled points, then
 * splines, normals, LE/TE, chord, TE parameters and panel angles. lconf: print a note. */
XFOIL_API void abcopy(xfoil_t *xf, int lconf);

#endif /* XFOIL_C_XGDES_H */
