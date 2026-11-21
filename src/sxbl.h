/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of s_xbl.f90 (routines split out of xbl to avoid circular module use).
 */
#ifndef XFOIL_C_SXBL_H
#define XFOIL_C_SXBL_H

#include "xfoil_c.h"
#include "xstate.h"

/* Sets the BL Newton-system line number isys for each BL station (from nbl)
 * and the system size nsys. Stops if nsys > 2*IVX. */
XFOIL_API void iblsys(xfoil_t *xf);

#endif /* XFOIL_C_SXBL_H */
