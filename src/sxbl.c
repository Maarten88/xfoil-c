/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of s_xbl.f90.
 */
#include "sxbl.h"

#include "xstate.h"
#include "xerror.h"

void iblsys(xfoil_t *xf)
{
    /* Fortran: do is = 1, 2; do ibl = 2, NBL(is); ISYS(ibl,is) = iv.
     * isys holds 0-based system line numbers (Fortran value - 1). */
    int iv = 0;         /* number of lines so far */
    for (int is = 0; is < 2; is++) {
        for (int ibl = 1; ibl < xf->nbl[is]; ibl++) {
            xf->isys[is][ibl] = iv;
            iv = iv + 1;
        }
    }

    xf->nsys = iv;
    if (xf->nsys > 2 * IVX) xfoil_stop("*** IBLSYS: BL system array overflow. ***");
}
