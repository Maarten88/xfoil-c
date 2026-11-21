/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of the API-reachable part of m_xfoil.f90.
 */
#ifndef XFOIL_C_XFOIL_H
#define XFOIL_C_XFOIL_H

#include "xfoil_c.h"
#include "xstate.h"

/* init (API entry point) is declared in xfoil_c.h. */

/* NACA 4- or 5-digit buffer airfoil ides1, then paneling (pangen). ides1 <= 0: Fortran
 * prompts for the designation (not ported): prints a message and returns. */
XFOIL_API void naca(xfoil_t *xf, int ides1);

/* Paneling of the current airfoil (x, y, n = npan) from the buffer airfoil (xb, yb, nb),
 * then TE, normals and panel angles. shopar: print the paneling parameters. */
XFOIL_API void pangen(xfoil_t *xf, int shopar);

#endif /* XFOIL_C_XFOIL_H */
