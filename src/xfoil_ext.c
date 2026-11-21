/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Extensions of the API (Phase 10b of MIGRATION_PLAN.md): functions that the original
 *  api.f90 does not have. They only call ported XFoil routines and do not change the
 *  state, so the drop-in API behaves exactly as before.
 */
#include <math.h>
#include <string.h>

#include "xstate.h"
#include "xerror.h"
#include "xfoil_c.h"
#include "xgeom.h"

static void get_geometry_impl(xfoil_t *xf, double *thick, double *x_thick, double *camber, double *x_camber, double *le_radius,
                  double *te_angle, double *area)
{
    if (xf->n < 2) {
        *thick = *x_thick = *camber = *x_camber = *le_radius = *te_angle = *area = NAN;
        return;
    }
    /* geopar as XFoil calls it for the buffer airfoil (set_airfoil), here on the current
     * airfoil; it writes only its arguments (t: work space). tccalc for the x positions,
     * which geopar computes but does not return. */
    double *t = xf->scr.geo_t;
    double sle, chord, ei11a, ei22a, apx1a, apx2a, ei11t, ei22t, apx1t, apx2t;
    geopar(xf->x, xf->xp, xf->y, xf->yp, xf->s, xf->n, t, &sle, &chord, area, le_radius, te_angle, &ei11a, &ei22a, &apx1a, &apx2a,
           &ei11t, &ei22t, &apx1t, &apx2t, thick, camber);
    double thick2, camber2;
    tccalc(xf->x, xf->xp, xf->y, xf->yp, xf->s, xf->n, &thick2, x_thick, &camber2, x_camber);
}

/* ======================================================================================
 * Flap deflection and hinge moment (Phase 10b item 2). Ported from the upstream XFOIL 6.99
 * source (xgdes.f FLAP, INSIDE, GETXYF; xgeom.f SSS, SCHECK), which xfoil-python does not
 * contain; plotting and console prompts dropped. Same exactness rules as the rest of the
 * port; the routines keep Fortran 1-based index values in their locals (IT1, IB1, ...).
 * ====================================================================================== */

#include <stdio.h>

#include "xfoil.h"
#include "xfortran.h"
#include "xspline.h"
#include "sxfoil.h"
#include "xblsys.h"

#define P2(v) __builtin_powi((v), 2)

/* INSIDE: true if (xf, yf) lies inside the closed contour x, y (n points). */
static int inside(const double *x_, const double *y_, int n_, double xf, double yf)
{
    /* integrate subtended angle around airfoil perimeter */
    double angle = 0.0;
    for (int i = 0; i < n_; i++) {
        const int ip = (i == n_ - 1) ? 0 : i + 1;
        const double xb1 = x_[i] - xf;
        const double yb1 = y_[i] - yf;
        const double xb2 = x_[ip] - xf;
        const double yb2 = y_[ip] - yf;
        angle = angle + (xb1 * yb2 - yb1 * xb2) / sqrt((P2(xb1) + P2(yb1)) * (P2(xb2) + P2(yb2)));
    }
    /* angle = 0 if xf,yf is outside, angle = +/- 2 pi if xf,yf is inside */
    return fabs(angle) > 1.0;
}

/* SSS: arc-length points s1 (fixed part) and s2 (flap) of the flap surface break; see xgeom.f. */
static void sss(double ss, double *s1, double *s2, double del, double xbf_, double ybf_, const double *x_,
                const double *xp_, const double *y_, const double *yp_, const double *s_, int n_, int iside)
{
    const double eps = 1.0e-5;      /* convergence epsilon */

    const double stot = fabs(s_[n_ - 1] - s_[0]);
    const double sind = sin(0.5 * fabs(del));

    double ssgn = 1.0;
    if (iside == 1) ssgn = -1.0;

    /* initial guesses for S1, S2 */
    const double rsq = P2(seval(ss, x_, xp_, s_, n_) - xbf_) + P2(seval(ss, y_, yp_, s_, n_) - ybf_);
    *s1 = ss - (sind * sqrt(rsq) + eps * stot) * ssgn;
    *s2 = ss + (sind * sqrt(rsq) + eps * stot) * ssgn;

    /* Newton iteration loop */
    int converged = 0;
    for (int iter = 1; iter <= 10; iter++) {
        const double x1 = seval(*s1, x_, xp_, s_, n_);
        const double x1p = deval(*s1, x_, xp_, s_, n_);
        const double y1 = seval(*s1, y_, yp_, s_, n_);
        const double y1p = deval(*s1, y_, yp_, s_, n_);

        const double x2 = seval(*s2, x_, xp_, s_, n_);
        const double x2p = deval(*s2, x_, xp_, s_, n_);
        const double y2 = seval(*s2, y_, yp_, s_, n_);
        const double y2p = deval(*s2, y_, yp_, s_, n_);

        const double r1sq = P2(x1 - xbf_) + P2(y1 - ybf_);
        const double r2sq = P2(x2 - xbf_) + P2(y2 - ybf_);
        const double r1 = sqrt(r1sq);
        const double r2 = sqrt(r2sq);

        const double rrsq = P2(x1 - x2) + P2(y1 - y2);
        const double rr = sqrt(rrsq);

        if (r1 <= eps * stot || r2 <= eps * stot) {
            *s1 = ss;
            *s2 = ss;
            return;
        }

        const double r1_s1 = (x1p * (x1 - xbf_) + y1p * (y1 - ybf_)) / r1;
        const double r2_s2 = (x2p * (x2 - xbf_) + y2p * (y2 - ybf_)) / r2;

        double rs1, a11, a12, rs2, a21, a22;
        if (sind > 0.01) {
            if (rr == 0.0) return;

            const double rr_s1 = (x1p * (x1 - x2) + y1p * (y1 - y2)) / rr;
            const double rr_s2 = -((x2p * (x1 - x2) + y2p * (y1 - y2)) / rr);

            /* residual 1: set included angle via dot product */
            rs1 = ((xbf_ - x1) * (x2 - x1) + (ybf_ - y1) * (y2 - y1)) / rr - sind * r1;
            a11 = ((xbf_ - x1) * (-x1p) + (ybf_ - y1) * (-y1p)) / rr
                + ((-x1p) * (x2 - x1) + (-y1p) * (y2 - y1)) / rr
                - ((xbf_ - x1) * (x2 - x1) + (ybf_ - y1) * (y2 - y1)) * rr_s1 / rrsq
                - sind * r1_s1;
            a12 = ((xbf_ - x1) * (x2p) + (ybf_ - y1) * (y2p)) / rr
                - ((xbf_ - x1) * (x2 - x1) + (ybf_ - y1) * (y2 - y1)) * rr_s2 / rrsq;

            /* residual 2: set equal length segments */
            rs2 = r1 - r2;
            a21 = r1_s1;
            a22 = -r2_s2;
        } else {
            /* residual 1: set included angle via small angle approximation */
            rs1 = (r1 + r2) * sind + (*s1 - *s2) * ssgn;
            a11 = r1_s1 * sind + ssgn;
            a12 = r2_s2 * sind - ssgn;

            /* residual 2: vector sum of the segments to the hinge perpendicular to the surface */
            const double x1pp = d2val(*s1, x_, xp_, s_, n_);
            const double y1pp = d2val(*s1, y_, yp_, s_, n_);
            const double x2pp = d2val(*s2, x_, xp_, s_, n_);
            const double y2pp = d2val(*s2, y_, yp_, s_, n_);

            const double xtot = x1 + x2 - 2.0 * xbf_;
            const double ytot = y1 + y2 - 2.0 * ybf_;

            rs2 = xtot * (x1p + x2p) + ytot * (y1p + y2p);
            a21 = x1p * (x1p + x2p) + y1p * (y1p + y2p) + xtot * x1pp + ytot * y1pp;
            a22 = x2p * (x1p + x2p) + y2p * (y1p + y2p) + xtot * x2pp + ytot * y2pp;
        }

        const double det = a11 * a22 - a12 * a21;
        double ds1 = -((rs1 * a22 - a12 * rs2) / det);
        double ds2 = -((a11 * rs2 - rs1 * a21) / det);

        ds1 = fmin_f(ds1, 0.01 * stot);
        ds1 = fmax_f(ds1, -.01 * stot);
        ds2 = fmin_f(ds2, 0.01 * stot);
        ds2 = fmax_f(ds2, -.01 * stot);

        *s1 = *s1 + ds1;
        *s2 = *s2 + ds2;
        if (fabs(ds1) + fabs(ds2) < eps * stot) {
            converged = 1;
            break;
        }
    }
    if (!converged) {
        if (xfoil_print) printf(" SSS: failed to converge subtending angle points\n");
        *s1 = ss;
        *s2 = ss;
    }

    /* make sure points are identical if included angle is zero */
    if (del == 0.0) {
        *s1 = 0.5 * (*s1 + *s2);
        *s2 = *s1;
    }
}

/* SCHECK: replaces the two nodes of a segment shorter than stol times an adjacent one by their
 * midpoint (xgeom.f). n may change. */
static void scheck(double *x_, double *y_, int *n_, double stol, int *lchange)
{
    *lchange = 0;
    /* check STOL for sanity */
    if (stol > 0.3) {
        if (xfoil_print) printf(" SCHECK:  Bad value for small panels (STOL > 0.3)\n");
        return;
    }
    /* restart the scan after every removal, until no segment is removed */
    int removed;
    do {
        removed = 0;
        for (int i = 1; i < *n_ - 2; i++) {
            const double dxm1 = x_[i] - x_[i - 1];
            const double dym1 = y_[i] - y_[i - 1];
            const double dsm1 = sqrt(dxm1 * dxm1 + dym1 * dym1);

            const double dxp1 = x_[i + 1] - x_[i];
            const double dyp1 = y_[i + 1] - y_[i];
            const double dsp1 = sqrt(dxp1 * dxp1 + dyp1 * dyp1);

            const double dxp2 = x_[i + 2] - x_[i + 1];
            const double dyp2 = y_[i + 2] - y_[i + 1];
            const double dsp2 = sqrt(dxp2 * dxp2 + dyp2 * dyp2);

            /* don't mess with doubled points (slope breaks) */
            if (dsp1 == 0.0) continue;

            if (dsp1 < stol * dsm1 || dsp1 < stol * dsp2) {
                /* replace node i with average of i and i+1 */
                x_[i] = 0.5 * (x_[i] + x_[i + 1]);
                y_[i] = 0.5 * (y_[i] + y_[i + 1]);
                /* remove node i+1 */
                for (int l = i + 1; l < *n_ - 1; l++) {       /* Fortran also read X(N+1) */
                    x_[l] = x_[l + 1];
                    y_[l] = y_[l + 1];
                }
                *n_ = *n_ - 1;
                *lchange = 1;
                if (xfoil_print) printf(" SCHECK segment removed at %12d\n", i + 1);
                removed = 1;
                break;
            }
        }
    } while (removed);
}

/* GETXYF without prompts: arc lengths tops/bots of the top/bottom surface at x = xf, and the
 * surface y values there (y/t = 1 and 0). */
static void getxyf(const double *x_, const double *xp_, const double *y_, const double *yp_, const double *s_,
                   int n_, double *tops, double *bots, double xf, double *topy, double *boty)
{
    /* find top and bottom y at hinge x location */
    *tops = s_[0] + (x_[0] - xf);
    *bots = s_[n_ - 1] - (x_[n_ - 1] - xf);
    sinvrt(tops, xf, x_, xp_, s_, n_);
    sinvrt(bots, xf, x_, xp_, s_, n_);
    *topy = seval(*tops, y_, yp_, s_, n_);
    *boty = seval(*bots, y_, yp_, s_, n_);

    if (xfoil_print) {
        printf("\n  Top    surface:  y =%8.4f     y/t = 1.0\n", *topy);
        printf("  Bottom surface:  y =%8.4f     y/t = 0.0\n", *boty);
    }
}


/* FLAP (xgdes.f) with RINPUT = (XBF, YBF, ddef), NINPUT = 3: deflects the buffer airfoil about
 * the hinge (xbf, ybf) by ddef degrees (+ = TE down). Points may be added/removed near the
 * surface breaks. ddef = 0: no change. */
static void flap_buffer(xfoil_t *xf, double ddef)
{
    double tops, bots, topy, boty;
    getxyf(xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb, &tops, &bots, xf->xbf, &topy, &boty);
    const int insid = inside(xf->xb, xf->yb, xf->nb, xf->xbf, xf->ybf);

    if (xfoil_print) printf("\n Flap hinge: x,y =%9.5f%9.5f\n", xf->xbf, xf->ybf);

    const double rdef = ddef * xf->pi / 180.0;
    if (rdef == 0.0) return;

    double atop, abot;
    if (insid) {
        atop = fmax_f(0.0, -rdef);
        abot = fmax_f(0.0, rdef);
    } else {
        const double chx = deval(bots, xf->xb, xf->xbp, xf->sb, xf->nb) - deval(tops, xf->xb, xf->xbp, xf->sb, xf->nb);
        const double chy = deval(bots, xf->yb, xf->ybp, xf->sb, xf->nb) - deval(tops, xf->yb, xf->ybp, xf->sb, xf->nb);
        const double fvx = seval(bots, xf->xb, xf->xbp, xf->sb, xf->nb) + seval(tops, xf->xb, xf->xbp, xf->sb, xf->nb);
        const double fvy = seval(bots, xf->yb, xf->ybp, xf->sb, xf->nb) + seval(tops, xf->yb, xf->ybp, xf->sb, xf->nb);
        const double crsp = chx * (xf->ybf - 0.5 * fvy) - chy * (xf->xbf - 0.5 * fvx);
        if (crsp > 0.0) {
            /* flap hinge is above airfoil */
            atop = fmax_f(0.0, rdef);
            abot = fmax_f(0.0, rdef);
        } else {
            /* flap hinge is below airfoil */
            atop = fmax_f(0.0, -rdef);
            abot = fmax_f(0.0, -rdef);
        }
    }

    /* find upper and lower surface break arc length values... */
    double st1, st2, sb1, sb2;
    sss(tops, &st1, &st2, atop, xf->xbf, xf->ybf, xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb, 1);
    sss(bots, &sb1, &sb2, abot, xf->xbf, xf->ybf, xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb, 2);

    /* ... and x,y coordinates */
    const double xt1 = seval(st1, xf->xb, xf->xbp, xf->sb, xf->nb);
    const double yt1 = seval(st1, xf->yb, xf->ybp, xf->sb, xf->nb);
    const double xt2 = seval(st2, xf->xb, xf->xbp, xf->sb, xf->nb);
    const double yt2 = seval(st2, xf->yb, xf->ybp, xf->sb, xf->nb);
    const double xb1 = seval(sb1, xf->xb, xf->xbp, xf->sb, xf->nb);
    const double yb1 = seval(sb1, xf->yb, xf->ybp, xf->sb, xf->nb);
    const double xb2 = seval(sb2, xf->xb, xf->xbp, xf->sb, xf->nb);
    const double yb2 = seval(sb2, xf->yb, xf->ybp, xf->sb, xf->nb);

    if (xfoil_print) {
        printf("\n Top breaks: x,y =  %9.5f%9.5f    %9.5f%9.5f\n", xt1, yt1, xt2, yt2);
        printf(" Bot breaks: x,y =  %9.5f%9.5f    %9.5f%9.5f\n", xb1, yb1, xb2, yb2);
    }

    /* find points adjacent to breaks (-1: not found) */
    int it1 = -1, it2 = -1, ib1 = -1, ib2 = -1;
    for (int i = 0; i < xf->nb - 1; i++) {
        if (xf->sb[i] <= st1 && xf->sb[i + 1] > st1) it1 = i + 1;
        if (xf->sb[i] < st2 && xf->sb[i + 1] >= st2) it2 = i;
        if (xf->sb[i] <= sb1 && xf->sb[i + 1] > sb1) ib1 = i;
        if (xf->sb[i] < sb2 && xf->sb[i + 1] >= sb2) ib2 = i + 1;
    }
    /* the Fortran left them undefined and then read outside the buffer (clean-up C6a) */
    if (it1 < 0 || it2 < 0 || ib1 < 1 || ib2 < 0 || it1 > xf->nb - 2)
        xfoil_stop("FLAP: surface break points not found (hinge outside the airfoil?)");

    const double dsavg = (xf->sb[xf->nb - 1] - xf->sb[0]) / (double)(xf->nb - 1);

    /* smallest fraction of s increments i+1 and i+2 away from break point */
    const double sfrac = 0.33333;

    double xt1new = 0.0, yt1new = 0.0, xt2new = 0.0, yt2new = 0.0;
    double xb1new = 0.0, yb1new = 0.0, xb2new = 0.0, yb2new = 0.0;
    int lt1new = 0, lt2new = 0, lb1new = 0, lb2new = 0;
    if (atop != 0.0) {
        const double st1p = st1 + sfrac * (xf->sb[it1] - st1);
        const double st1q = st1 + sfrac * (xf->sb[it1 + 1] - st1);
        if (xf->sb[it1] < st1q) {
            /* simply move adjacent point to ideal SFRAC location */
            xt1new = seval(st1q, xf->xb, xf->xbp, xf->sb, xf->nb);
            yt1new = seval(st1q, xf->yb, xf->ybp, xf->sb, xf->nb);
            lt1new = 0;
        } else {
            /* make new point at SFRAC location */
            xt1new = seval(st1p, xf->xb, xf->xbp, xf->sb, xf->nb);
            yt1new = seval(st1p, xf->yb, xf->ybp, xf->sb, xf->nb);
            lt1new = 1;
        }

        const double st2p = st2 + sfrac * (xf->sb[it2] - st2);
        const int it2q = (it2 > 1) ? it2 - 1 : 0;
        const double st2q = st2 + sfrac * (xf->sb[it2q] - st2);
        if (xf->sb[it2] > st2q) {
            /* simply move adjacent point */
            xt2new = seval(st2q, xf->xb, xf->xbp, xf->sb, xf->nb);
            yt2new = seval(st2q, xf->yb, xf->ybp, xf->sb, xf->nb);
            lt2new = 0;
        } else {
            /* make new point */
            xt2new = seval(st2p, xf->xb, xf->xbp, xf->sb, xf->nb);
            yt2new = seval(st2p, xf->yb, xf->ybp, xf->sb, xf->nb);
            lt2new = 1;
        }
    }

    if (abot != 0.0) {
        const double sb1p = sb1 + sfrac * (xf->sb[ib1] - sb1);
        const double sb1q = sb1 + sfrac * (xf->sb[ib1 - 1] - sb1);
        if (xf->sb[ib1] > sb1q) {
            /* simply move adjacent point */
            xb1new = seval(sb1q, xf->xb, xf->xbp, xf->sb, xf->nb);
            yb1new = seval(sb1q, xf->yb, xf->ybp, xf->sb, xf->nb);
            lb1new = 0;
        } else {
            /* make new point */
            xb1new = seval(sb1p, xf->xb, xf->xbp, xf->sb, xf->nb);
            yb1new = seval(sb1p, xf->yb, xf->ybp, xf->sb, xf->nb);
            lb1new = 1;
        }

        const double sb2p = sb2 + sfrac * (xf->sb[ib2] - sb2);
        const int ib2q = (ib2 + 1 < xf->nb - 1) ? ib2 + 1 : xf->nb - 1;
        const double sb2q = sb2 + sfrac * (xf->sb[ib2q] - sb2);
        if (xf->sb[ib2] < sb2q) {
            /* simply move adjacent point */
            xb2new = seval(sb2q, xf->xb, xf->xbp, xf->sb, xf->nb);
            yb2new = seval(sb2q, xf->yb, xf->ybp, xf->sb, xf->nb);
            lb2new = 0;
        } else {
            /* make new point */
            xb2new = seval(sb2p, xf->xb, xf->xbp, xf->sb, xf->nb);
            yb2new = seval(sb2p, xf->yb, xf->ybp, xf->sb, xf->nb);
            lb2new = 1;
        }
    }

    const double sind = sin(rdef);
    const double cosd = cos(rdef);

    /* rotate flap points about the hinge point (XBF,YBF) */
    for (int i = 0; i < xf->nb; i++) {
        if (i >= it1 && i <= ib1) continue;

        const double xbar = xf->xb[i] - xf->xbf;
        const double ybar = xf->yb[i] - xf->ybf;

        xf->xb[i] = xf->xbf + xbar * cosd + ybar * sind;
        xf->yb[i] = xf->ybf - xbar * sind + ybar * cosd;
    }

    int idif = it1 - it2 - 1;
    if (idif > 0) {
        /* delete points on upper airfoil surface which "disappeared" */
        xf->nb = xf->nb - idif;
        it1 = it1 - idif;
        ib1 = ib1 - idif;
        ib2 = ib2 - idif;
        for (int i = it2 + 1; i < xf->nb; i++) {
            xf->sb[i] = xf->sb[i + idif];
            xf->xb[i] = xf->xb[i + idif];
            xf->yb[i] = xf->yb[i + idif];
        }
    }

    idif = ib2 - ib1 - 1;
    if (idif > 0) {
        /* delete points on lower airfoil surface which "disappeared" */
        xf->nb = xf->nb - idif;
        ib2 = ib2 - idif;
        for (int i = ib1 + 1; i < xf->nb; i++) {
            xf->sb[i] = xf->sb[i + idif];
            xf->xb[i] = xf->xb[i + idif];
            xf->yb[i] = xf->yb[i + idif];
        }
    }

    if (atop == 0.0) {
        /* arc length of newly created surface on top of airfoil */
        const double dsnew = fabs(rdef) * sqrt(P2(xt1 - xf->xbf) + P2(yt1 - xf->ybf));

        /* number of points to be added to define newly created surface */
        const int npadd = (int)(1.5 * dsnew / dsavg + 1.0);

        /* skip everything if no points are to be added */
        if (npadd != 0) {
            /* increase coordinate array length to make room for the new point(s) */
            xf->nb = xf->nb + npadd;
            it1 = it1 + npadd;
            ib1 = ib1 + npadd;
            ib2 = ib2 + npadd;
            for (int i = xf->nb - 1; i >= it1; i--) {
                xf->xb[i] = xf->xb[i - npadd];
                xf->yb[i] = xf->yb[i - npadd];
            }

            /* add new points along the new surface circular arc segment */
            const double dang = rdef / (double)npadd;
            const double xbar = xt1 - xf->xbf;
            const double ybar = yt1 - xf->ybf;
            for (int ip = 1; ip <= npadd; ip++) {
                const double ang = dang * ((double)ip - 0.5);
                const double ca = cos(ang);
                const double sa = sin(ang);

                xf->xb[it1 - ip] = xf->xbf + xbar * ca + ybar * sa;
                xf->yb[it1 - ip] = xf->ybf - xbar * sa + ybar * ca;
            }
        }
    } else {
        /* set point in the corner and possibly two adjacent points */
        int npadd = 1;
        if (lt2new) npadd = npadd + 1;
        if (lt1new) npadd = npadd + 1;

        xf->nb = xf->nb + npadd;
        it1 = it1 + npadd;
        ib1 = ib1 + npadd;
        ib2 = ib2 + npadd;
        for (int i = xf->nb - 1; i >= it1; i--) {
            xf->xb[i] = xf->xb[i - npadd];
            xf->yb[i] = xf->yb[i - npadd];
        }

        if (lt1new) {
            xf->xb[it1 - 1] = xt1new;
            xf->yb[it1 - 1] = yt1new;
            xf->xb[it1 - 2] = xt1;
            xf->yb[it1 - 2] = yt1;
        } else {
            xf->xb[it1] = xt1new;
            xf->yb[it1] = yt1new;
            xf->xb[it1 - 1] = xt1;
            xf->yb[it1 - 1] = yt1;
        }

        const double xbar = xt2new - xf->xbf;
        const double ybar = yt2new - xf->ybf;
        if (lt2new) {
            xf->xb[it2 + 1] = xf->xbf + xbar * cosd + ybar * sind;
            xf->yb[it2 + 1] = xf->ybf - xbar * sind + ybar * cosd;
        } else {
            xf->xb[it2] = xf->xbf + xbar * cosd + ybar * sind;
            xf->yb[it2] = xf->ybf - xbar * sind + ybar * cosd;
        }
    }

    if (abot == 0.0) {
        /* arc length of newly created surface on bottom of airfoil */
        const double dsnew = fabs(rdef) * sqrt(P2(xb1 - xf->xbf) + P2(yb1 - xf->ybf));

        /* number of points to be added to define newly created surface */
        const int npadd = (int)(1.5 * dsnew / dsavg + 1.0);

        /* skip everything if no points are to be added */
        if (npadd != 0) {
            /* increase coordinate array length to make room for the new point(s) */
            xf->nb = xf->nb + npadd;
            ib2 = ib2 + npadd;
            for (int i = xf->nb - 1; i >= ib2; i--) {
                xf->xb[i] = xf->xb[i - npadd];
                xf->yb[i] = xf->yb[i - npadd];
            }

            /* add new points along the new surface circular arc segment */
            const double dang = rdef / (double)npadd;
            const double xbar = xb1 - xf->xbf;
            const double ybar = yb1 - xf->ybf;
            for (int ip = 1; ip <= npadd; ip++) {
                const double ang = dang * ((double)ip - 0.5);
                const double ca = cos(ang);
                const double sa = sin(ang);

                xf->xb[ib1 + ip] = xf->xbf + xbar * ca + ybar * sa;
                xf->yb[ib1 + ip] = xf->ybf - xbar * sa + ybar * ca;
            }
        }
    } else {
        /* set point in the corner and possibly two adjacent points */
        int npadd = 1;
        if (lb2new) npadd = npadd + 1;
        if (lb1new) npadd = npadd + 1;

        xf->nb = xf->nb + npadd;
        ib2 = ib2 + npadd;
        for (int i = xf->nb - 1; i >= ib2; i--) {
            xf->xb[i] = xf->xb[i - npadd];
            xf->yb[i] = xf->yb[i - npadd];
        }

        if (lb1new) {
            xf->xb[ib1 + 1] = xb1new;
            xf->yb[ib1 + 1] = yb1new;
            xf->xb[ib1 + 2] = xb1;
            xf->yb[ib1 + 2] = yb1;
        } else {
            xf->xb[ib1] = xb1new;
            xf->yb[ib1] = yb1new;
            xf->xb[ib1 + 1] = xb1;
            xf->yb[ib1 + 1] = yb1;
        }

        const double xbar = xb2new - xf->xbf;
        const double ybar = yb2new - xf->ybf;
        if (lb2new) {
            xf->xb[ib2 - 1] = xf->xbf + xbar * cosd + ybar * sind;
            xf->yb[ib2 - 1] = xf->ybf - xbar * sind + ybar * cosd;
        } else {
            xf->xb[ib2] = xf->xbf + xbar * cosd + ybar * sind;
            xf->yb[ib2] = xf->ybf - xbar * sind + ybar * cosd;
        }
    }


    /* check new geometry for splinter segments */
    const double stol = 0.2;
    int lchange;
    scheck(xf->xb, xf->yb, &xf->nb, stol, &lchange);

    /* spline new geometry */
    scalc(xf->xb, xf->yb, xf->sb, xf->nb);
    segspl(xf->xb, xf->xbp, xf->sb, xf->nb);
    segspl(xf->yb, xf->ybp, xf->sb, xf->nb);

    geopar(xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb, xf->w1, &xf->sble, &xf->chordb, &xf->areab, &xf->radble, &xf->angbte, &xf->ei11ba, &xf->ei22ba, &xf->apx1ba,
           &xf->apx2ba, &xf->ei11bt, &xf->ei22bt, &xf->apx1bt, &xf->apx2bt, &xf->thickb, &xf->cambrb);

    xf->lbflap = 1;

    if (xf->lgsym) {
        if (xfoil_print) {
            printf("\n");
            printf(" Disabling symmetry enforcement\n");
        }
        xf->lgsym = 0;
    }
    /* Fortran: plotting of the hinge, breaks and new shape follows (dropped) */
}


/* ---- API ------------------------------------------------------------------------------ */

/* The undeflected buffer airfoil, kept at the first flap after a new airfoil. */
static void flap_save_base(xfoil_t *xf)
{
    memcpy(xf->scr.flap_base_xb, xf->xb, sizeof xf->xb);
    memcpy(xf->scr.flap_base_yb, xf->yb, sizeof xf->yb);
    xf->scr.flap_base_nb = xf->nb;
    xf->scr.flap_base_valid = 1;
}

/* Buffer airfoil := undeflected airfoil, with its splines and geometry parameters (as FLAP leaves them). */
static void flap_restore_base(xfoil_t *xf)
{
    memcpy(xf->xb, xf->scr.flap_base_xb, sizeof xf->xb);
    memcpy(xf->yb, xf->scr.flap_base_yb, sizeof xf->yb);
    xf->nb = xf->scr.flap_base_nb;
    scalc(xf->xb, xf->yb, xf->sb, xf->nb);
    segspl(xf->xb, xf->xbp, xf->sb, xf->nb);
    segspl(xf->yb, xf->ybp, xf->sb, xf->nb);
    geopar(xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb, xf->w1, &xf->sble, &xf->chordb, &xf->areab, &xf->radble,
           &xf->angbte, &xf->ei11ba, &xf->ei22ba, &xf->apx1ba, &xf->apx2ba, &xf->ei11bt, &xf->ei22bt, &xf->apx1bt,
           &xf->apx2bt, &xf->thickb, &xf->cambrb);
}

/* Undeflected airfoil deflected by `angle` (absolute) about the current hinge, then repaneled. */
static void flap_apply(xfoil_t *xf, double angle)
{
    if (xf->scr.flap_base_valid)
        flap_restore_base(xf);
    else
        flap_save_base(xf);
    xf->scr.flap_angle = angle;
    if (angle != 0.0) flap_buffer(xf, angle);
    pangen(xf, 1);                          /* repanel with the current parameters; copies the hinge */
    xf->hmom = xf->hfx = xf->hfy = NAN;
}

static void set_hinge_impl(xfoil_t *xf, const double *x_hinge, const double *y_hinge, const bool *y_relative)
{
    if (xf->nb < 2) return;
    /* the hinge is defined on the undeflected airfoil */
    if (xf->scr.flap_base_valid) flap_restore_base(xf);
    double tops, bots, topy, boty;
    getxyf(xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb, &tops, &bots, *x_hinge, &topy, &boty);
    const double yf = *y_relative ? topy * *y_hinge + boty * (1.0 - *y_hinge) : *y_hinge;

    /* buffer hinge (FLAP, copied to the current airfoil by pangen) and current hinge (FMOM) */
    xf->xbf = *x_hinge;
    xf->ybf = yf;
    xf->lbflap = 1;
    xf->xof = xf->xbf;
    xf->yof = xf->ybf;
    xf->lflap = 1;
    xf->hmom = xf->hfx = xf->hfy = NAN;     /* until the next analysis */

    /* a deflected flap moves with the hinge: re-apply the current angle about the new hinge */
    if (xf->scr.flap_base_valid && xf->scr.flap_angle != 0.0) {
        xf->scr.cp_valid = 0;
        flap_apply(xf, xf->scr.flap_angle);
    }
}

static void flap_impl(xfoil_t *xf, const double *angle)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */
    if (!xf->lbflap || xf->nb < 2) return;
    /* absolute deflection (user decision 2026-11-21): the undeflected airfoil deflected by angle;
     * XFOIL's FLAP command adds to the current deflection */
    flap_apply(xf, *angle);
}

double get_flap_angle(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    return xf->scr.flap_angle;
}

int get_hinge_point(xfoil_t *xf, double *x_hinge, double *y_hinge)
{
    xfoil_print = xf->show_output;
    *x_hinge = xf->lflap ? xf->xof : NAN;
    *y_hinge = xf->lflap ? xf->yof : NAN;
    return xf->lflap != 0;
}

void get_hinge(xfoil_t *xf, double *hinge_moment, double *hinge_fx, double *hinge_fy)
{
    xfoil_print = xf->show_output;
    *hinge_moment = xf->lflap ? xf->hmom : NAN;
    *hinge_fx = xf->lflap ? xf->hfx : NAN;
    *hinge_fy = xf->lflap ? xf->hfy : NAN;
}

void clear_hinge(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    xf->xbf = 0.0;
    xf->ybf = 0.0;
    xf->lbflap = 0;
    xf->xof = 0.0;
    xf->yof = 0.0;
    xf->lflap = 0;
    xf->hmom = xf->hfx = xf->hfy = NAN;
    xf->scr.flap_base_valid = 0;            /* a new airfoil: no deflection to undo */
    xf->scr.flap_angle = 0.0;
}

/* ---- public entry points: run inside the error guard (xerror.h) ---- */

void get_geometry(xfoil_t *xf, double *thick, double *x_thick, double *camber, double *x_camber, double *le_radius,
                  double *te_angle, double *area)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(get_geometry_impl(xf, thick, x_thick, camber, x_camber, le_radius, te_angle, area));
}

void set_hinge(xfoil_t *xf, const double *x_hinge, const double *y_hinge, const bool *y_relative)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(set_hinge_impl(xf, x_hinge, y_hinge, y_relative));
}

void flap(xfoil_t *xf, const double *angle)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(flap_impl(xf, angle));
}

/* ---- boundary layer of the last analysed point (XFOIL OPER DUMP, extension) ---- */

int get_n_bl(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    return xf->n + (xf->lwake ? xf->nw : 0);
}

/* The quantities of XFOIL's bldump for the airfoil nodes (i < n, in airfoil order) and the wake nodes
 * (n <= i < n + nw), with the same formulas; NaN (except s, x, y) unless the last analysed point is a
 * converged viscous one. ampl: amplification factor N on laminar stations, ctau: max shear stress
 * coefficient on turbulent stations and in the wake (the other one NaN; XFOIL keeps both in CTAU).
 * Wake: cf = 0 and hstar NaN, as bldump writes no H* there. */
static void get_bl_impl(xfoil_t *xf, double *s, double *x, double *y, double *ue, double *dstar, double *theta,
                        double *cf, double *hk, double *hstar, double *ampl, double *ctau_out, const int *n_points)
{
    const int nn = xf->n + (xf->lwake ? xf->nw : 0);
    const int n_out = *n_points < nn ? *n_points : nn;
    const int valid = xf->scr.cp_valid && xf->lipan && xf->lvisc;
    comset(xf);
    const double hstinv = xf->gamm1 * P2(xf->minf / xf->qinf) / (1.0 + 0.5 * xf->gamm1 * P2(xf->minf));

    for (int i = 0; i < n_out; i++) {
        s[i] = xf->s[i];
        x[i] = xf->x[i];
        y[i] = xf->y[i];
        ue[i] = dstar[i] = theta[i] = cf[i] = hk[i] = hstar[i] = ampl[i] = ctau_out[i] = NAN;
        if (!valid) continue;

        const int wake = i >= xf->n;
        const int is = (wake || xf->gam[i] < 0.0) ? 1 : 0;
        const int ibl = is == 0 ? xf->iblte[0] - i : xf->iblte[1] + i + 1 - xf->n;
        const double ds = xf->dstr[is][ibl];
        const double th = xf->thet[is][ibl];
        double h, hs;
        if (wake) {
            h = ds / th;
            hs = NAN;
            cf[i] = 0.0;
        } else if (th == 0.0) {
            h = 1.0;
            hs = 1.0;
            cf[i] = xf->tau[is][ibl] / (0.5 * P2(xf->qinf));
        } else {
            h = ds / th;
            hs = xf->tstr[is][ibl] / th;
            cf[i] = xf->tau[is][ibl] / (0.5 * P2(xf->qinf));
        }
        const double g = wake ? xf->uedg[is][ibl] : xf->gam[i];
        const double u = (g / xf->qinf) * (1.0 - xf->tklam) / (1.0 - xf->tklam * P2(g / xf->qinf));
        const double amsq = u * u * hstinv / (xf->gamm1 * (1.0 - 0.5 * u * u * hstinv));
        double dummy;
        hkin(h, amsq, &hk[i], &dummy, &dummy);
        ue[i] = u;
        dstar[i] = ds;
        theta[i] = th;
        hstar[i] = hs;
        if (ibl < xf->itran[is]) ampl[i] = xf->ctau[is][ibl];
        else ctau_out[i] = xf->ctau[is][ibl];
    }
}

void get_bl(xfoil_t *xf, double *s, double *x, double *y, double *ue, double *dstar, double *theta, double *cf,
            double *hk, double *hstar, double *ampl, double *ctau, const int *n_points)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(get_bl_impl(xf, s, x, y, ue, dstar, theta, cf, hk, hstar, ampl, ctau, n_points));
}
