/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of s_xoper.f90.
 */
#include "sxoper.h"

#include "xstate.h"
#include "xspline.h"

/* Fortran S(i), X(i), ... with 1-based i (literal translation of the panel loops). */

void mhinge(xfoil_t *xf)
{
    double tops, bots;
    const double *cp = xf->lvisc ? xf->cpv : xf->cpi;

    if (!xf->lflap) {
        /* Fortran: call getxyf(...) (console prompt for the hinge), LFLap = .true.
         * Not ported (no console prompts); unreachable, callers check lflap. */
        return;
    }
    /* find top and bottom y at hinge x location */
    tops = xf->xof;
    bots = xf->s[xf->n - 1] - xf->xof;
    sinvrt(&tops, xf->xof, xf->x, xf->xp, xf->s, xf->n);
    sinvrt(&bots, xf->xof, xf->x, xf->xp, xf->s, xf->n);

    const double topx = seval(tops, xf->x, xf->xp, xf->s, xf->n);
    const double topy = seval(tops, xf->y, xf->yp, xf->s, xf->n);
    const double botx = seval(bots, xf->x, xf->xp, xf->s, xf->n);
    const double boty = seval(bots, xf->y, xf->yp, xf->s, xf->n);

    xf->hmom = 0.0;
    xf->hfx = 0.0;
    xf->hfy = 0.0;

    double dx, dy, xmid, ymid, pmid, frac;
    int i;

    /* integrate pressures on top and bottom sides of flap */
    for (i = 1; i < xf->n; i++) {
        if (xf->s[i - 1] < tops || xf->s[i] > bots) {
            dx = xf->x[i] - xf->x[i - 1];
            dy = xf->y[i] - xf->y[i - 1];
            xmid = 0.5 * (xf->x[i] + xf->x[i - 1]) - xf->xof;
            ymid = 0.5 * (xf->y[i] + xf->y[i - 1]) - xf->yof;
            pmid = 0.5 * (cp[i] + cp[i - 1]);
            xf->hmom = xf->hmom + pmid * (xmid * dx + ymid * dy);
            xf->hfx = xf->hfx - pmid * dy;
            xf->hfy = xf->hfy + pmid * dx;
        }
    }

    /* find the s[i-1]..s[i] interval containing s=TOPS (no hit: last panel; the Fortran read
     * S(N+1) beyond the airfoil there; clean-up C6a) */
    for (i = 1; i < xf->n; i++) {
        if (xf->s[i] > tops) break;
    }
    if (i > xf->n - 1) i = xf->n - 1;

    /* add on top surface chunk TOPS..s[i-1], missed in the loop above */
    dx = topx - xf->x[i - 1];
    dy = topy - xf->y[i - 1];
    xmid = 0.5 * (topx + xf->x[i - 1]) - xf->xof;
    ymid = 0.5 * (topy + xf->y[i - 1]) - xf->yof;
    if (xf->s[i] != xf->s[i - 1]) {
        frac = (tops - xf->s[i - 1]) / (xf->s[i] - xf->s[i - 1]);
    } else {
        frac = 0.0;
    }
    const double topp = cp[i] * frac + cp[i - 1] * (1.0 - frac);
    pmid = 0.5 * (topp + cp[i - 1]);
    xf->hmom = xf->hmom + pmid * (xmid * dx + ymid * dy);
    xf->hfx = xf->hfx - pmid * dy;
    xf->hfy = xf->hfy + pmid * dx;

    /* add on inside flap surface contribution from hinge to top surface */
    dx = xf->xof - topx;
    dy = xf->yof - topy;
    xmid = 0.5 * (topx + xf->xof) - xf->xof;
    ymid = 0.5 * (topy + xf->yof) - xf->yof;
    xf->hmom = xf->hmom + pmid * (xmid * dx + ymid * dy);
    xf->hfx = xf->hfx - pmid * dy;
    xf->hfy = xf->hfy + pmid * dx;

    /* find the s[i-1]..s[i] interval containing s=BOTS. Fortran ends with I = 1 (C 0)
     * when there is no hit and then reads S(0) (out of bounds); C stops at i = 1 instead
     * (Phase 12). */
    for (i = xf->n - 1; i >= 1; i--) {
        if (xf->s[i - 1] < bots) break;
    }
    if (i < 1) i = 1;

    /* add on bottom surface chunk BOTS..s[i], missed in the loop above */
    dx = xf->x[i] - botx;
    dy = xf->y[i] - boty;
    xmid = 0.5 * (botx + xf->x[i]) - xf->xof;
    ymid = 0.5 * (boty + xf->y[i]) - xf->yof;
    if (xf->s[i] != xf->s[i - 1]) {
        frac = (bots - xf->s[i - 1]) / (xf->s[i] - xf->s[i - 1]);
    } else {
        frac = 0.0;
    }
    const double botp = cp[i] * frac + cp[i - 1] * (1.0 - frac);
    pmid = 0.5 * (botp + cp[i]);
    xf->hmom = xf->hmom + pmid * (xmid * dx + ymid * dy);
    xf->hfx = xf->hfx - pmid * dy;
    xf->hfy = xf->hfy + pmid * dx;

    /* add on inside flap surface contribution from hinge to bottom surface */
    dx = botx - xf->xof;
    dy = boty - xf->yof;
    xmid = 0.5 * (botx + xf->xof) - xf->xof;
    ymid = 0.5 * (boty + xf->yof) - xf->yof;
    xf->hmom = xf->hmom + pmid * (xmid * dx + ymid * dy);
    xf->hfx = xf->hfx - pmid * dy;
    xf->hfy = xf->hfy + pmid * dx;

    /* add on TE base thickness contribution */
    dx = xf->x[0] - xf->x[xf->n - 1];
    dy = xf->y[0] - xf->y[xf->n - 1];
    xmid = 0.5 * (xf->x[0] + xf->x[xf->n - 1]) - xf->xof;
    ymid = 0.5 * (xf->y[0] + xf->y[xf->n - 1]) - xf->yof;
    pmid = 0.5 * (cp[0] + cp[xf->n - 1]);
    xf->hmom = xf->hmom + pmid * (xmid * dx + ymid * dy);
    xf->hfx = xf->hfx - pmid * dy;
    xf->hfy = xf->hfy + pmid * dx;
}
