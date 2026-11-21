/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of m_xutils.f90. Fortran `real**integer` is evaluated by gfortran with
 *  __builtin_powi (repeated squaring), not pow(); it is used here too so that
 *  results are bit-identical.
 */
#include "xutils.h"

#include <math.h>
#include <stdio.h>

#include "xstate.h"
#include "xerror.h"
#include "xfortran.h"

void setexp(double *s, double ds1, double smax, int nn)
{
    const double sigma = smax / ds1;
    const int nex = nn - 1;
    const double rnex = (double)nex;
    const double rni = 1.0 / rnex;

    /* solve quadratic for initial geometric ratio guess */
    const double aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
    const double bbb = rnex * (rnex - 1.0) / 2.0;
    const double ccc = rnex - sigma;

    double disc = bbb * bbb - 4.0 * aaa * ccc;
    disc = fmax_f(0.0, disc);

    double ratio;
    if (nex <= 1) {
        xfoil_stop("SETEXP: Cannot fill array.  N too small.");
    } else if (nex == 2) {
        ratio = -(ccc / bbb) + 1.0;
    } else {
        ratio = (-bbb + sqrt(disc)) / (2.0 * aaa) + 1.0;
    }

    if (ratio != 1.0) {
        /* Newton iteration for actual geometric ratio */
        int converged = 0;
        for (int iter = 1; iter <= 100; iter++) {
            const double sigman = (__builtin_powi(ratio, nex) - 1.0) / (ratio - 1.0);
            const double res = pow(sigman, rni) - pow(sigma, rni);
            const double dresdr = rni * pow(sigman, rni)
                                  * (rnex * __builtin_powi(ratio, nex - 1) - sigman)
                                  / (__builtin_powi(ratio, nex) - 1.0);

            const double dratio = -(res / dresdr);
            ratio = ratio + dratio;

            if (fabs(dratio) < 1.0e-5) {
                converged = 1;
                break;
            }
        }
        if (!converged && xfoil_print)
            printf(" SETEXP: Convergence failed.  Continuing anyway ...\n");
    }

    /* set up stretched array using converged geometric ratio */
    s[0] = 0.0;
    double ds = ds1;
    for (int n = 1; n < nn; n++) {
        s[n] = s[n - 1] + ds;
        ds = ds * ratio;
    }
}

double atanc(double y, double x, double thold)
{
    static const double pi = 3.1415926535897932384;
    static const double tpi = 6.2831853071795864769;

    /* new position angle, ignoring the atan2 branch cut for now */
    const double thnew = atan2(y, x);
    const double dthet = thnew - thold;

    /* the angle change cannot exceed +/- pi: remove multiples of 2 pi.
     * Fortran int() truncates to a default (32-bit) integer; sign() keeps -0.0. */
    const double dtcorr = dthet - tpi * (double)(int)((dthet + copysign(pi, dthet)) / tpi);

    return thold + dtcorr;
}
