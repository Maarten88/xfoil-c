/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of m_spline.f90. Expression and evaluation order follow the Fortran
 *  (e.g. `a*b/c` is (a*b)/c, `x**2` is x*x), so results are bit-identical.
 */
#include "xspline.h"

#include <math.h>
#include <stdio.h>

#include "xstate.h"
#include "xerror.h"
#include "xfortran.h"

#define NMAX 1000   /* scratch size of splind, as in the Fortran */

/* Binary search used by seval/deval/d2val/curv: returns i (0-based, 1..n-1)
 * with s[i-1] <= ss < s[i] (clamped to the end intervals). Same steps as the
 * Fortran 1-based search: every imid is the Fortran imid minus one. */
static int find_interval(double ss, const double *s, int n)
{
    int ilow = 0;
    int i = n - 1;
    while (i - ilow > 1) {
        const int imid = (i + ilow) / 2;
        if (ss < s[imid])
            i = imid;
        else
            ilow = imid;
    }
    return i;
}

void splind(const double *x, double *xs, const double *s, int n, double xs1, double xs2)
{
    double a[NMAX], b[NMAX], c[NMAX];

    if (n > NMAX) xfoil_stop("SPLIND: array overflow, increase NMAX");

    for (int i = 1; i < n - 1; i++) {
        const double dsm = s[i] - s[i - 1];
        const double dsp = s[i + 1] - s[i];
        b[i] = dsp;
        a[i] = 2.0 * (dsm + dsp);
        c[i] = dsm;
        xs[i] = 3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
    }

    if (xs1 == 999.0) {
        /* zero second derivative */
        a[0] = 2.0;
        c[0] = 1.0;
        xs[0] = 3.0 * (x[1] - x[0]) / (s[1] - s[0]);
    } else if (xs1 == -999.0) {
        /* zero third derivative */
        a[0] = 1.0;
        c[0] = 1.0;
        xs[0] = 2.0 * (x[1] - x[0]) / (s[1] - s[0]);
    } else {
        /* specified first derivative */
        a[0] = 1.0;
        c[0] = 0.0;
        xs[0] = xs1;
    }

    if (xs2 == 999.0) {
        b[n - 1] = 1.0;
        a[n - 1] = 2.0;
        xs[n - 1] = 3.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);
    } else if (xs2 == -999.0) {
        b[n - 1] = 1.0;
        a[n - 1] = 1.0;
        xs[n - 1] = 2.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);
    } else {
        a[n - 1] = 1.0;
        b[n - 1] = 0.0;
        xs[n - 1] = xs2;
    }

    if (n == 2 && xs1 == -999.0 && xs2 == -999.0) {
        b[n - 1] = 1.0;
        a[n - 1] = 2.0;
        xs[n - 1] = 3.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);
    }

    trisol(a, b, c, xs, n);
}

void trisol(double *a, const double *b, double *c, double *d, int kk)
{
    for (int k = 1; k < kk; k++) {
        const int km = k - 1;
        c[km] = c[km] / a[km];
        d[km] = d[km] / a[km];
        a[k] = a[k] - b[k] * c[km];
        d[k] = d[k] - b[k] * d[km];
    }

    d[kk - 1] = d[kk - 1] / a[kk - 1];

    for (int k = kk - 2; k >= 0; k--)
        d[k] = d[k] - c[k] * d[k + 1];
}

double seval(double ss, const double *x, const double *xs, const double *s, int n)
{
    const int i = find_interval(ss, s, n);
    const double ds = s[i] - s[i - 1];
    const double t = (ss - s[i - 1]) / ds;
    const double cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    const double cx2 = ds * xs[i] - x[i] + x[i - 1];
    return t * x[i] + (1.0 - t) * x[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2);
}

double deval(double ss, const double *x, const double *xs, const double *s, int n)
{
    const int i = find_interval(ss, s, n);
    const double ds = s[i] - s[i - 1];
    const double t = (ss - s[i - 1]) / ds;
    const double cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    const double cx2 = ds * xs[i] - x[i] + x[i - 1];
    double d = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    return d / ds;
}

double d2val(double ss, const double *x, const double *xs, const double *s, int n)
{
    const int i = find_interval(ss, s, n);
    const double ds = s[i] - s[i - 1];
    const double t = (ss - s[i - 1]) / ds;
    const double cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    const double cx2 = ds * xs[i] - x[i] + x[i - 1];
    double d = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    return d / (ds * ds);
}

double curv(double ss, const double *x, const double *xs, const double *y,
            const double *ys, const double *s, int n)
{
    const int i = find_interval(ss, s, n);
    const double ds = s[i] - s[i - 1];
    const double t = (ss - s[i - 1]) / ds;

    const double cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    const double cx2 = ds * xs[i] - x[i] + x[i - 1];
    const double xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    const double xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;

    const double cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
    const double cy2 = ds * ys[i] - y[i] + y[i - 1];
    const double yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
    const double ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;

    double sd = sqrt(xd * xd + yd * yd);
    const double sdmin = 0.001 * ds;
    sd = fmax_f(sd, sdmin);             /* Fortran max(sd, 0.001*ds) */

    return (xd * ydd - yd * xdd) / (sd * sd * sd);
}

void sinvrt(double *si, double xi, const double *x, const double *xs, const double *s, int n)
{
    const double sisav = *si;

    for (int iter = 1; iter <= 10; iter++) {
        const double res = seval(*si, x, xs, s, n) - xi;
        const double resp = deval(*si, x, xs, s, n);
        const double ds = -(res / resp);
        *si = *si + ds;
        if (fabs(ds / (s[n - 1] - s[0])) < 1.0e-5) return;
    }
    if (xfoil_print) printf(" SINVRT: spline inversion failed. Input value returned.\n");
    *si = sisav;
}

void scalc(const double *x, const double *y, double *s, int n)
{
    s[0] = 0.0;
    for (int i = 1; i < n; i++) {
        const double dx = x[i] - x[i - 1];
        const double dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + sqrt(dx * dx + dy * dy);
    }
}

void segspl(const double *x, double *xs, const double *s, int n)
{
    if (s[0] == s[1]) xfoil_stop("SEGSPL:  First input point duplicated");
    if (s[n - 1] == s[n - 2]) xfoil_stop("SEGSPL:  Last  input point duplicated");

    int iseg0 = 0;
    for (int iseg = 1; iseg <= n - 3; iseg++) {
        if (s[iseg] == s[iseg + 1]) {
            const int nseg = iseg - iseg0 + 1;
            splind(x + iseg0, xs + iseg0, s + iseg0, nseg, -999.0, -999.0);
            iseg0 = iseg + 1;
        }
    }

    const int nseg = n - iseg0;
    splind(x + iseg0, xs + iseg0, s + iseg0, nseg, -999.0, -999.0);
}

