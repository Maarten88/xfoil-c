/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  NACA 4- and 5-digit airfoil generators, ported from m_naca.f90.
 *  Conventions as in xsort.h; strings as in xuserio.h (blank-padded, explicit length).
 *  Exported (XFOIL_API) so the tests can compare them with the Fortran originals.
 */
#ifndef XFOIL_C_XNACA_H
#define XFOIL_C_XNACA_H

#include "xfoil_c.h"

/* NACA 4-digit airfoil `ides` (0..9999) with nside points per side.
 * Work arrays xx, yt, yc: nside each. Output: xb, yb (2*nside), *nb = 2*nside-1
 * points from the upper-surface TE round the LE to the lower-surface TE, and
 * name = "NACA dddd" blank-padded to name_len (>= 9). */
XFOIL_API void naca4(int ides, double *xx, double *yt, double *yc, int nside,
                     double *xb, double *yb, int *nb, char *name, int name_len);

/* NACA 5-digit airfoil; the first three digits must be 210..250 in steps of 10.
 * Otherwise prints a message (show_output), sets *ides = 0 and returns without
 * touching the outputs. name_len >= 10. */
XFOIL_API void naca5(int *ides, double *xx, double *yt, double *yc, int nside,
                     double *xb, double *yb, int *nb, char *name, int name_len);

#endif /* XFOIL_C_XNACA_H */
