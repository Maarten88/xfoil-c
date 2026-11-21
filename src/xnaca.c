/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of m_naca.f90. Integer powers use __builtin_powi, as gfortran does for
 *  x**n (also for constant n: x**4 is (x*x)*(x*x)); real powers use pow().
 */
#include "xnaca.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "xstate.h"

static const double AN = 1.5;           /* TE point bunching parameter */
static const char DIGITS[] = "0123456789";

/* x/c distribution with TE bunching (shared by naca4 and naca5) */
static double xx_point(int i, int nside)    /* i is 1-based, as in the Fortran */
{
    const double anp = AN + 1.0;
    const double frac = (double)(i - 1) / (double)(nside - 1);
    if (i == nside) return 1.0;
    return 1.0 - anp * frac * pow(1.0 - frac, AN) - pow(1.0 - frac, anp);
}

/* NACA 4-digit thickness distribution */
static double thickness(double x, double t)
{
    return (0.29690 * sqrt(x) - 0.12600 * x - 0.35160 * __builtin_powi(x, 2)
            + 0.28430 * __builtin_powi(x, 3) - 0.10150 * __builtin_powi(x, 4)) * t / 0.20;
}

/* Upper surface TE -> LE, then lower surface LE -> TE */
static void assemble(const double *xx, const double *yt, const double *yc, int nside,
                     double *xb, double *yb, int *nb)
{
    int ib = 0;
    for (int i = nside - 1; i >= 0; i--) {
        xb[ib] = xx[i];
        yb[ib] = yc[i] + yt[i];
        ib++;
    }
    for (int i = 1; i < nside; i++) {
        xb[ib] = xx[i];
        yb[ib] = yc[i] - yt[i];
        ib++;
    }
    *nb = ib;
}

/* Fortran: Name = 'NACA'; Name(6:5+ndig) = digits */
static void set_name(char *name, int name_len, const int *dig, int ndig)
{
    memset(name, ' ', (size_t)name_len);
    memcpy(name, "NACA", 4);
    for (int k = 0; k < ndig; k++) name[5 + k] = DIGITS[dig[k]];
}

void naca4(int ides, double *xx, double *yt, double *yc, int nside,
           double *xb, double *yb, int *nb, char *name, int name_len)
{
    const int n4 = ides / 1000;
    const int n3 = (ides - n4 * 1000) / 100;
    const int n2 = (ides - n4 * 1000 - n3 * 100) / 10;
    const int n1 = (ides - n4 * 1000 - n3 * 100 - n2 * 10);

    const double m = (double)n4 / 100.0;
    const double p = (double)n3 / 10.0;
    const double t = (double)(n2 * 10 + n1) / 100.0;

    for (int i = 0; i < nside; i++) {
        const double x = xx[i] = xx_point(i + 1, nside);
        yt[i] = thickness(x, t);
        if (x < p)
            yc[i] = m / __builtin_powi(p, 2) * (2.0 * p * x - __builtin_powi(x, 2));
        else
            yc[i] = m / __builtin_powi(1.0 - p, 2) * ((1.0 - 2.0 * p) + 2.0 * p * x - __builtin_powi(x, 2));
    }

    assemble(xx, yt, yc, nside, xb, yb, nb);

    const int dig[4] = {n4, n3, n2, n1};
    set_name(name, name_len, dig, 4);
}

void naca5(int *ides, double *xx, double *yt, double *yc, int nside,
           double *xb, double *yb, int *nb, char *name, int name_len)
{
    const int d = *ides;
    const int n5 = d / 10000;
    const int n4 = (d - n5 * 10000) / 1000;
    const int n3 = (d - n5 * 10000 - n4 * 1000) / 100;
    const int n2 = (d - n5 * 10000 - n4 * 1000 - n3 * 100) / 10;
    const int n1 = (d - n5 * 10000 - n4 * 1000 - n3 * 100 - n2 * 10);

    const int n543 = 100 * n5 + 10 * n4 + n3;

    double m, c;
    switch (n543) {
    case 210: m = 0.0580; c = 361.4;  break;    /* P = 0.05 */
    case 220: m = 0.1260; c = 51.64;  break;    /* P = 0.10 */
    case 230: m = 0.2025; c = 15.957; break;    /* P = 0.15 */
    case 240: m = 0.2900; c = 6.643;  break;    /* P = 0.20 */
    case 250: m = 0.3910; c = 3.230;  break;    /* P = 0.25 */
    default:
        if (xfoil_print) {
            printf(" Illegal 5-digit designation\n");
            printf(" First three digits must be 210, 220, ... 250\n");
        }
        *ides = 0;
        return;
    }

    const double t = (double)(n2 * 10 + n1) / 100.0;

    for (int i = 0; i < nside; i++) {
        const double x = xx[i] = xx_point(i + 1, nside);
        yt[i] = thickness(x, t);
        if (x < m)
            yc[i] = (c / 6.0) * (__builtin_powi(x, 3) - 3.0 * m * __builtin_powi(x, 2) + m * m * (3.0 - m) * x);
        else
            yc[i] = (c / 6.0) * __builtin_powi(m, 3) * (1.0 - x);
    }

    assemble(xx, yt, yc, nside, xb, yb, nb);

    const int dig[5] = {n5, n4, n3, n2, n1};
    set_name(name, name_len, dig, 5);
}
