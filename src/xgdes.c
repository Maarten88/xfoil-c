/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of the API-reachable part of m_xgdes.f90.
 */
#include "xgdes.h"

#include <math.h>
#include <stdio.h>

#include "xstate.h"
#include "sxfoil.h"
#include "xgeom.h"
#include "xpanel.h"
#include "xspline.h"

#define P2(v) __builtin_powi((v), 2)
/* Fortran 1-based accessors (literal translation of the doubled-point loop) */

void abcopy(xfoil_t *xf, int lconf)
{
    if (xf->nb <= 1) {
        if (xfoil_print) printf(" ABCOPY: Buffer airfoil not available.\n");
        return;
    } else if (xf->nb > IQX - 5) {
        if (xfoil_print) {
            printf(" Maximum number of panel nodes  : %12d\n", IQX - 5);
            printf(" Number of buffer airfoil points: %12d\n", xf->nb);
            printf(" Current airfoil cannot be set.\n");
            printf(" Try executing PANE at Top Level instead.\n");
        }
        return;
    }
    if (xf->n != xf->nb) xf->lblini = 0;

    xf->n = xf->nb;
    for (int i = 0; i < xf->n; i++) {
        xf->x[i] = xf->xb[i];
        xf->y[i] = xf->yb[i];
    }

    if (xf->lbflap) {
        xf->xof = xf->xbf;
        xf->yof = xf->ybf;
        xf->lflap = 1;
    }

    /* strip out doubled points (the index is not re-checked after a removal, so a
     * tripled point leaves a doubled one, as in Fortran) */
    int i = 0;
    for (;;) {
        i = i + 1;
        if (xf->x[i - 1] == xf->x[i] && xf->y[i - 1] == xf->y[i]) {
            for (int j = i; j < xf->n - 1; j++) {
                xf->x[j] = xf->x[j + 1];
                xf->y[j] = xf->y[j + 1];
            }
            xf->n = xf->n - 1;
        }
        if (i >= xf->n - 1) break;
    }

    scalc(xf->x, xf->y, xf->s, xf->n);
    segspl(xf->x, xf->xp, xf->s, xf->n);
    segspl(xf->y, xf->yp, xf->s, xf->n);

    ncalc(xf->x, xf->y, xf->s, xf->n, xf->nx, xf->ny);

    lefind(&xf->sle, xf->x, xf->xp, xf->y, xf->yp, xf->s, xf->n);
    xf->xle = seval(xf->sle, xf->x, xf->xp, xf->s, xf->n);
    xf->yle = seval(xf->sle, xf->y, xf->yp, xf->s, xf->n);
    xf->xte = 0.5 * (xf->x[0] + xf->x[xf->n - 1]);
    xf->yte = 0.5 * (xf->y[0] + xf->y[xf->n - 1]);
    xf->chord = sqrt(P2(xf->xte - xf->xle) + P2(xf->yte - xf->yle));

    tecalc(xf);
    apcalc(xf);

    xf->lgamu = 0;
    xf->lwake = 0;
    xf->lqaij = 0;
    xf->ladij = 0;
    xf->lwdij = 0;
    xf->lipan = 0;
    xf->lvconv = 0;

    if (lconf && xfoil_print) printf("\n Current airfoil nodes set from buffer airfoil nodes (%4d )\n", xf->n);
}
