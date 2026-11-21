/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of the API-reachable routines of m_xgeom.f90.
 */
#include "xgeom.h"

#include <math.h>
#include <stdio.h>

#include "xstate.h"
#include "xfortran.h"
#include "xspline.h"
#include "xutils.h"

void lefind(double *sle, const double *x, const double *xp, const double *y,
            const double *yp, const double *s, int n)
{
    /* convergence tolerance */
    const double dseps = (s[n - 1] - s[0]) * 1.0e-5;

    /* trailing edge point */
    const double xte = 0.5 * (x[0] + x[n - 1]);
    const double yte = 0.5 * (y[0] + y[n - 1]);

    /* first guess: first point where the surface turns back towards the TE.
     * Fortran do i = 3, N-2: after a full loop i = N-1, i.e. 0-based n-2;
     * after a zero-trip loop i = 3, i.e. 0-based 2. Same for this for-loop. */
    int i;
    for (i = 2; i <= n - 3; i++) {
        const double dxte = x[i] - xte;
        const double dyte = y[i] - yte;
        const double dx = x[i + 1] - x[i];
        const double dy = y[i + 1] - y[i];
        const double dotp = dxte * dx + dyte * dy;
        if (dotp < 0.0) break;
    }

    *sle = s[i];

    /* sharp LE */
    if (s[i] == s[i - 1]) return;

    /* Newton iteration for the exact SLE */
    for (int iter = 1; iter <= 50; iter++) {
        const double xle = seval(*sle, x, xp, s, n);
        const double yle = seval(*sle, y, yp, s, n);
        const double dxds = deval(*sle, x, xp, s, n);
        const double dyds = deval(*sle, y, yp, s, n);
        const double dxdd = d2val(*sle, x, xp, s, n);
        const double dydd = d2val(*sle, y, yp, s, n);

        const double xchord = xle - xte;
        const double ychord = yle - yte;

        /* drive dot product between chord line and LE tangent to zero */
        const double res = xchord * dxds + ychord * dyds;
        const double ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;

        double dsle = -(res / ress);
        dsle = fmax_f(dsle, -(0.02 * fabs(xchord + ychord)));
        dsle = fmin_f(dsle, 0.02 * fabs(xchord + ychord));
        *sle = *sle + dsle;
        if (fabs(dsle) < dseps) return;
    }
    if (xfoil_print) printf(" LEFIND:  LE point not found.  Continuing...\n");
    *sle = s[i];
}

void sopps(double *sopp, double si, const double *x, const double *xp, const double *y,
           const double *yp, const double *s, int n, double sle)
{
    /* reference length for testing convergence */
    const double slen = s[n - 1] - s[0];

    /* chord line unit vector */
    double xle = seval(sle, x, xp, s, n);
    double yle = seval(sle, y, yp, s, n);
    const double xte = 0.5 * (x[0] + x[n - 1]);
    const double yte = 0.5 * (y[0] + y[n - 1]);
    const double chord = sqrt(__builtin_powi(xte - xle, 2) + __builtin_powi(yte - yle, 2));
    const double dxc = (xte - xle) / chord;
    const double dyc = (yte - yle) / chord;

    int in, inopp;
    if (si < sle) {
        in = 0;
        inopp = n - 1;
    } else {
        in = n - 1;
        inopp = 0;
    }
    const double sfrac = (si - sle) / (s[in] - sle);
    *sopp = sle + sfrac * (s[inopp] - sle);

    if (fabs(sfrac) <= 1.0e-5) {
        *sopp = sle;
        return;
    }

    /* xbar = x coordinate in chord-line axes */
    const double xi = seval(si, x, xp, s, n);
    const double yi = seval(si, y, yp, s, n);
    xle = seval(sle, x, xp, s, n);
    yle = seval(sle, y, yp, s, n);
    const double xbar = (xi - xle) * dxc + (yi - yle) * dyc;

    /* converge on the opposite point with the same xbar */
    for (int iter = 1; iter <= 12; iter++) {
        const double xopp = seval(*sopp, x, xp, s, n);
        const double yopp = seval(*sopp, y, yp, s, n);
        const double xoppd = deval(*sopp, x, xp, s, n);
        const double yoppd = deval(*sopp, y, yp, s, n);

        const double res = (xopp - xle) * dxc + (yopp - yle) * dyc - xbar;
        const double resd = xoppd * dxc + yoppd * dyc;

        if (fabs(res) / slen < 1.0e-5) return;
        if (resd == 0.0) break;

        const double dsopp = -(res / resd);
        *sopp = *sopp + dsopp;

        if (fabs(dsopp) / slen < 1.0e-5) return;
    }
    if (xfoil_print) printf(" SOPPS: Opposite-point location failed. Continuing...\n");
    *sopp = sle + sfrac * (s[inopp] - sle);
}

void norm(double *x, double *xp, double *y, double *yp, double *s, int n)
{
    double sle;

    scalc(x, y, s, n);
    segspl(x, xp, s, n);
    segspl(y, yp, s, n);

    lefind(&sle, x, xp, y, yp, s, n);

    const double xmax = 0.5 * (x[0] + x[n - 1]);
    const double xmin = seval(sle, x, xp, s, n);
    const double ymin = seval(sle, y, yp, s, n);

    const double fudge = 1.0 / (xmax - xmin);
    for (int i = 0; i < n; i++) {
        x[i] = (x[i] - xmin) * fudge;
        y[i] = (y[i] - ymin) * fudge;
        s[i] = s[i] * fudge;
    }
}

void geopar(const double *x, const double *xp, const double *y, const double *yp,
            const double *s, int n, double *t, double *sle, double *chord, double *area,
            double *radle, double *angte, double *ei11a, double *ei22a, double *apx1a,
            double *apx2a, double *ei11t, double *ei22t, double *apx1t, double *apx2t,
            double *thick, double *cambr)
{
    double xcena, ycena, xcent, ycent, slen, xthick, xcambr;

    lefind(sle, x, xp, y, yp, s, n);

    const double xle = seval(*sle, x, xp, s, n);
    const double yle = seval(*sle, y, yp, s, n);
    const double xte = 0.5 * (x[0] + x[n - 1]);
    const double yte = 0.5 * (y[0] + y[n - 1]);

    const double chsq = __builtin_powi(xte - xle, 2) + __builtin_powi(yte - yle, 2);
    *chord = sqrt(chsq);

    const double curvle = curv(*sle, x, xp, y, yp, s, n);

    *radle = 0.0;
    if (fabs(curvle) > 0.001 * (s[n - 1] - s[0])) *radle = 1.0 / curvle;

    const double ang1 = atan2(-yp[0], -xp[0]);
    const double ang2 = atanc(yp[n - 1], xp[n - 1], ang1);
    *angte = ang2 - ang1;

    for (int i = 0; i < n; i++) t[i] = 1.0;

    aecalc(n, x, y, t, 1, area, &xcena, &ycena, ei11a, ei22a, apx1a, apx2a);
    aecalc(n, x, y, t, 2, &slen, &xcent, &ycent, ei11t, ei22t, apx1t, apx2t);

    /* approximate thickness and camber (on the discrete points only) */
    tccalc(x, xp, y, yp, s, n, thick, &xthick, cambr, &xcambr);

    if (xfoil_print) {
        printf(" Max thickness = %12.6f  at x = %7.3f\n Max camber    = %12.6f  at x = %7.3f\n", *thick, xthick, *cambr, xcambr);
    }
}

void aecalc(int n, const double *x, const double *y, const double *t, int itype,
            double *area, double *xcen, double *ycen, double *ei11, double *ei22,
            double *apx1, double *apx2)
{
    static const double pi = 3.141592653589793238;

    double sint = 0.0, aint = 0.0, xint = 0.0, yint = 0.0;
    double xxint = 0.0, xyint = 0.0, yyint = 0.0;

    for (int io = 0; io < n; io++) {
        const int ip = (io == n - 1) ? 0 : io + 1;

        const double dx = x[io] - x[ip];
        const double dy = y[io] - y[ip];
        const double xa = (x[io] + x[ip]) * 0.50;
        const double ya = (y[io] + y[ip]) * 0.50;
        const double ta = (t[io] + t[ip]) * 0.50;

        const double ds = sqrt(dx * dx + dy * dy);
        sint = sint + ds;

        if (itype == 1) {
            /* integrate over the airfoil cross-section */
            const double da = ya * dx;
            aint = aint + da;
            xint = xint + xa * da;
            yint = yint + ya * da / 2.0;
            xxint = xxint + xa * xa * da;
            xyint = xyint + xa * ya * da / 2.0;
            yyint = yyint + ya * ya * da / 3.0;
        } else {
            /* integrate over the skin thickness */
            const double da = ta * ds;
            aint = aint + da;
            xint = xint + xa * da;
            yint = yint + ya * da;
            xxint = xxint + xa * xa * da;
            xyint = xyint + xa * ya * da;
            yyint = yyint + ya * ya * da;
        }
    }

    *area = aint;

    if (aint == 0.0) {
        *xcen = 0.0;
        *ycen = 0.0;
        *ei11 = 0.0;
        *ei22 = 0.0;
        *apx1 = 0.0;
        *apx2 = atan2(1.0, 0.0);
        return;
    }

    /* centroid */
    *xcen = xint / aint;
    *ycen = yint / aint;

    /* inertias */
    const double eixx = yyint - *ycen * *ycen * aint;
    const double eixy = xyint - *xcen * *ycen * aint;
    const double eiyy = xxint - *xcen * *xcen * aint;

    /* principal-axis inertias; ei11 is closest to "up-down" bending inertia */
    const double eisq = 0.25 * __builtin_powi(eixx - eiyy, 2) + __builtin_powi(eixy, 2);
    const double sgn = copysign(1.0, eiyy - eixx);
    *ei11 = 0.5 * (eixx + eiyy) - sgn * sqrt(eisq);
    *ei22 = 0.5 * (eixx + eiyy) + sgn * sqrt(eisq);

    if (*ei11 == 0.0 || *ei22 == 0.0) {
        /* vanishing section stiffness */
        *apx1 = 0.0;
        *apx2 = atan2(1.0, 0.0);
    } else if (eisq / (*ei11 * *ei22) < __builtin_powi(0.001 * sint, 4)) {
        /* rotationally-invariant section (circle, square, etc.) */
        *apx1 = 0.0;
        *apx2 = atan2(1.0, 0.0);
    } else {
        /* normal airfoil section */
        const double c1 = eixy;
        const double s1 = eixx - *ei11;
        const double c2 = eixy;
        const double s2 = eixx - *ei22;

        if (fabs(s1) > fabs(s2)) {
            *apx1 = atan2(s1, c1);
            *apx2 = *apx1 + 0.5 * pi;
        } else {
            *apx2 = atan2(s2, c2);
            *apx1 = *apx2 - 0.5 * pi;
        }

        if (*apx1 < -0.5 * pi) *apx1 = *apx1 + pi;
        if (*apx1 > +0.5 * pi) *apx1 = *apx1 - pi;
        if (*apx2 < -0.5 * pi) *apx2 = *apx2 + pi;
        if (*apx2 > +0.5 * pi) *apx2 = *apx2 - pi;
    }
}

void tccalc(const double *x, const double *xp, const double *y, const double *yp,
            const double *s, int n, double *thick, double *xthick, double *cambr,
            double *xcambr)
{
    double sle, sopp;

    lefind(&sle, x, xp, y, yp, s, n);
    const double xle = seval(sle, x, xp, s, n);
    const double yle = seval(sle, y, yp, s, n);
    const double xte = 0.5 * (x[0] + x[n - 1]);
    const double yte = 0.5 * (y[0] + y[n - 1]);
    const double chord = sqrt(__builtin_powi(xte - xle, 2) + __builtin_powi(yte - yle, 2));

    /* unit chord-line vector */
    const double dxc = (xte - xle) / chord;
    const double dyc = (yte - yle) / chord;

    *thick = 0.0;
    *xthick = 0.0;
    *cambr = 0.0;
    *xcambr = 0.0;

    /* y-thickness and camber at each point */
    for (int i = 0; i < n; i++) {
        const double ybar = (y[i] - yle) * dxc - (x[i] - xle) * dyc;

        /* point on the opposite side with the same chord x value */
        sopps(&sopp, s[i], x, xp, y, yp, s, n, sle);
        const double xopp = seval(sopp, x, xp, s, n);
        const double yopp = seval(sopp, y, yp, s, n);

        const double ybarop = (yopp - yle) * dxc - (xopp - xle) * dyc;

        const double yc = 0.5 * (ybar + ybarop);
        const double yt = fabs(ybar - ybarop);

        if (fabs(yc) > fabs(*cambr)) {
            *cambr = yc;
            *xcambr = xopp;
        }
        if (fabs(yt) > fabs(*thick)) {
            *thick = yt;
            *xthick = xopp;
        }
    }
}

void cang(const double *x, const double *y, int n, int iprint, int *imax, double *amax)
{
    *amax = 0.0;
    *imax = 0;

    if (iprint == 2 && xfoil_print) printf("\n  i       x        y      angle\n");

    for (int i = 1; i < n - 1; i++) {
        double dx1 = x[i] - x[i - 1];
        double dy1 = y[i] - y[i - 1];
        double dx2 = x[i] - x[i + 1];
        double dy2 = y[i] - y[i + 1];

        /* allow for doubled points; a doubled point at an end has no panel on that side, so
         * there is no corner (Fortran read X(0) / X(N+1) there; clean-up C6a) */
        if (dx1 == 0.0 && dy1 == 0.0) {
            if (i < 2) continue;
            dx1 = x[i] - x[i - 2];
            dy1 = y[i] - y[i - 2];
        }
        if (dx2 == 0.0 && dy2 == 0.0) {
            if (i + 2 > n - 1) continue;
            dx2 = x[i] - x[i + 2];
            dy2 = y[i] - y[i + 2];
        }

        const double crossp = (dx2 * dy1 - dy2 * dx1)
                              / sqrt((__builtin_powi(dx1, 2) + __builtin_powi(dy1, 2))
                                     * (__builtin_powi(dx2, 2) + __builtin_powi(dy2, 2)));
        const double angl = asin(crossp) * (180.0 / 3.1415926);

        if (iprint == 2 && xfoil_print) {
            printf(" %3d%9.4f%9.4f%9.3f\n", i + 1, x[i], y[i], angl);     /* index printed 1-based */
        }
        if (fabs(angl) > fabs(*amax)) {
            *amax = angl;
            *imax = i;
        }
    }

    if (iprint >= 1 && xfoil_print) {
        printf("\n Maximum panel corner angle =%7.3f   at  i,x,y  = %3d%9.4f%9.4f\n", *amax, *imax + 1, x[*imax], y[*imax]);
    }
}
