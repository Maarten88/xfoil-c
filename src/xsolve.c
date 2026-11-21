/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of ludcmp and baksub from m_xsolve.f90.
 */
#include "xsolve.h"

#include <math.h>
#include <stddef.h>

#include "xstate.h"
#include "xerror.h"
#include "xfortran.h"

#define NVX 500

/* Fortran A(i,j), 0-based, leading dimension nsiz */
#define A(i, j) a[(size_t)(j) * (size_t)nsiz + (size_t)(i)]

void ludcmp(int nsiz, int n, double *a, int *indx)
{
    double vv[NVX];
    int imax = 0;       /* uninitialised in the Fortran; only read if a pivot is NaN */

    if (n > NVX) xfoil_stop("LUDCMP: Array overflow. Increase NVX.");

    for (int i = 0; i < n; i++) {
        double aamax = 0.0;
        for (int j = 0; j < n; j++) aamax = fmax_f(fabs(A(i, j)), aamax);
        vv[i] = 1.0 / aamax;
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < j; i++) {
            double sum = A(i, j);
            for (int k = 0; k < i; k++) sum = sum - A(i, k) * A(k, j);
            A(i, j) = sum;
        }

        double aamax = 0.0;
        for (int i = j; i < n; i++) {
            double sum = A(i, j);
            for (int k = 0; k < j; k++) sum = sum - A(i, k) * A(k, j);
            A(i, j) = sum;

            const double dum = vv[i] * fabs(sum);
            if (dum >= aamax) {
                imax = i;
                aamax = dum;
            }
        }

        if (j != imax) {
            for (int k = 0; k < n; k++) {
                const double dum = A(imax, k);
                A(imax, k) = A(j, k);
                A(j, k) = dum;
            }
            vv[imax] = vv[j];
        }

        indx[j] = imax;
        if (j != n - 1) {
            const double dum = 1.0 / A(j, j);
            for (int i = j + 1; i < n; i++) A(i, j) = A(i, j) * dum;
        }
    }
}

void baksub(int nsiz, int n, const double *a, const int *indx, double *b)
{
    /* forward substitution; ii = first non-zero element of b (-1: none yet,
     * Fortran ii = 0) */
    int ii = -1;
    for (int i = 0; i < n; i++) {
        const int ll = indx[i];
        double sum = b[ll];
        b[ll] = b[i];
        if (ii != -1) {
            for (int j = ii; j < i; j++) sum = sum - A(i, j) * b[j];
        } else if (sum != 0.0) {
            ii = i;
        }
        b[i] = sum;
    }

    /* back substitution */
    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];
        if (i < n - 1) {
            for (int j = i + 1; j < n; j++) sum = sum - A(i, j) * b[j];
        }
        b[i] = sum / A(i, i);
    }
}

void gauss(int nsiz, int nn, double *z, double *r, int nrhs)
{
#define Z(i, j) z[(size_t)(j) * (size_t)nsiz + (size_t)(i)]
#define R(i, j) r[(size_t)(j) * (size_t)nsiz + (size_t)(i)]
    for (int np = 0; np < nn - 1; np++) {
        const int np1 = np + 1;

        /* max pivot index nx */
        int nx = np;
        for (int n = np1; n < nn; n++)
            if (fabs(Z(n, np)) > fabs(Z(nx, np))) nx = n;

        const double pivot = 1.0 / Z(nx, np);

        /* switch pivots */
        Z(nx, np) = Z(np, np);

        /* switch rows and normalize the pivot row */
        for (int l = np1; l < nn; l++) {
            const double temp = Z(nx, l) * pivot;
            Z(nx, l) = Z(np, l);
            Z(np, l) = temp;
        }
        for (int l = 0; l < nrhs; l++) {
            const double temp = R(nx, l) * pivot;
            R(nx, l) = R(np, l);
            R(np, l) = temp;
        }

        /* forward eliminate everything */
        for (int k = np1; k < nn; k++) {
            const double ztmp = Z(k, np);
            for (int l = np1; l < nn; l++) Z(k, l) = Z(k, l) - ztmp * Z(np, l);
            for (int l = 0; l < nrhs; l++) R(k, l) = R(k, l) - ztmp * R(np, l);
        }
    }

    /* solve for the last row */
    for (int l = 0; l < nrhs; l++) R(nn - 1, l) = R(nn - 1, l) / Z(nn - 1, nn - 1);

    /* back substitute everything */
    for (int np = nn - 2; np >= 0; np--) {
        const int np1 = np + 1;
        for (int l = 0; l < nrhs; l++)
            for (int k = np1; k < nn; k++) R(np, l) = R(np, l) - Z(np, k) * R(k, l);
    }
#undef Z
#undef R
}

#if XFOIL_FAST
/* Row updates of the fast blsolv over contiguous, non-overlapping rows (restrict: vectorisable).
 * Same arithmetic as the exact loops: y[l] - a*x[l], y[l] - (a*x[l] + b*z[l] [+ c*w[l]]). */
static inline void row_sub1(double *restrict y, const double *restrict x, double a, int lo, int hi)
{
    for (int l = lo; l < hi; l++) y[l] = y[l] - a * x[l];
}

static inline void row_sub2(double *restrict y, const double *restrict x, const double *restrict z,
                            double a, double b, int lo, int hi)
{
    for (int l = lo; l < hi; l++) y[l] = y[l] - (a * x[l] + b * z[l]);
}

static inline void row_sub3(double *restrict y, const double *restrict x, const double *restrict z,
                            const double *restrict w, double a, double b, double c, int lo, int hi)
{
    for (int l = lo; l < hi; l++) y[l] = y[l] - (a * x[l] + b * z[l] + c * w[l]);
}

void blsolv(xfoil_t *xf)
{
    /* XFOIL_FAST (optimization O2): VM rows are contiguous (xlayout.h), so every loop over l below
     * runs over contiguous memory and vectorises; the arithmetic is that of the exact version. */
    double (*restrict vm)[3][IZX] = (double (*)[3][IZX])(void *)xf->vm;
    const int nsys = xf->nsys;         /* local: writes through vm cannot change it (vectorisation) */
    /* Blocks, all 0-based: Fortran VA(k,l,iv) is va[iv][l-1][k-1] (row k, column l),
     * VM(k,jv,iv) is vm[iv][jv][k-1], VDEL(k,l,iv) is vdel[iv][l-1][k-1]. */

    const int ivte1 = xf->isys[0][xf->iblte[0]];

    const double vacc1 = xf->vaccel;
    const double vacc2 = xf->vaccel * 2.0 / (xf->s[xf->n - 1] - xf->s[0]);
    const double vacc3 = xf->vaccel * 2.0 / (xf->s[xf->n - 1] - xf->s[0]);

    for (int iv = 0; iv < nsys; iv++) {
        const int ivp = iv + 1;

        /* ==== invert the VA(iv) block */

        /* normalize the first row */
        double pivot = 1.0 / xf->va[iv][0][0];
        xf->va[iv][1][0] = xf->va[iv][1][0] * pivot;
        for (int l = iv; l < nsys; l++) vm[iv][0][l] = vm[iv][0][l] * pivot;
        xf->vdel[iv][0][0] = xf->vdel[iv][0][0] * pivot;
        xf->vdel[iv][1][0] = xf->vdel[iv][1][0] * pivot;

        /* eliminate the lower first column in the VA block */
        for (int k = 1; k < 3; k++) {
            const double vtmp = xf->va[iv][0][k];
            xf->va[iv][1][k] = xf->va[iv][1][k] - vtmp * xf->va[iv][1][0];
            for (int l = iv; l < nsys; l++) vm[iv][k][l] = vm[iv][k][l] - vtmp * vm[iv][0][l];
            xf->vdel[iv][0][k] = xf->vdel[iv][0][k] - vtmp * xf->vdel[iv][0][0];
            xf->vdel[iv][1][k] = xf->vdel[iv][1][k] - vtmp * xf->vdel[iv][1][0];
        }

        /* normalize the second row */
        pivot = 1.0 / xf->va[iv][1][1];
        for (int l = iv; l < nsys; l++) vm[iv][1][l] = vm[iv][1][l] * pivot;
        xf->vdel[iv][0][1] = xf->vdel[iv][0][1] * pivot;
        xf->vdel[iv][1][1] = xf->vdel[iv][1][1] * pivot;

        /* eliminate the lower second column in the VA block */
        {
            const int k = 2;   /* third row */
            const double vtmp = xf->va[iv][1][k];
            for (int l = iv; l < nsys; l++) vm[iv][k][l] = vm[iv][k][l] - vtmp * vm[iv][1][l];
            xf->vdel[iv][0][k] = xf->vdel[iv][0][k] - vtmp * xf->vdel[iv][0][1];
            xf->vdel[iv][1][k] = xf->vdel[iv][1][k] - vtmp * xf->vdel[iv][1][1];
        }

        /* normalize the third row */
        pivot = 1.0 / vm[iv][2][iv];
        for (int l = ivp; l < nsys; l++) vm[iv][2][l] = vm[iv][2][l] * pivot;
        xf->vdel[iv][0][2] = xf->vdel[iv][0][2] * pivot;
        xf->vdel[iv][1][2] = xf->vdel[iv][1][2] * pivot;

        /* eliminate the upper third column in the VA block */
        {
            const double vtmp1 = vm[iv][0][iv];
            const double vtmp2 = vm[iv][1][iv];
            for (int l = ivp; l < nsys; l++) {
                vm[iv][0][l] = vm[iv][0][l] - vtmp1 * vm[iv][2][l];
                vm[iv][1][l] = vm[iv][1][l] - vtmp2 * vm[iv][2][l];
            }
            xf->vdel[iv][0][0] = xf->vdel[iv][0][0] - vtmp1 * xf->vdel[iv][0][2];
            xf->vdel[iv][0][1] = xf->vdel[iv][0][1] - vtmp2 * xf->vdel[iv][0][2];
            xf->vdel[iv][1][0] = xf->vdel[iv][1][0] - vtmp1 * xf->vdel[iv][1][2];
            xf->vdel[iv][1][1] = xf->vdel[iv][1][1] - vtmp2 * xf->vdel[iv][1][2];
        }

        /* eliminate the upper second column in the VA block */
        {
            const double vtmp = xf->va[iv][1][0];
            for (int l = ivp; l < nsys; l++) vm[iv][0][l] = vm[iv][0][l] - vtmp * vm[iv][1][l];
            xf->vdel[iv][0][0] = xf->vdel[iv][0][0] - vtmp * xf->vdel[iv][0][1];
            xf->vdel[iv][1][0] = xf->vdel[iv][1][0] - vtmp * xf->vdel[iv][1][1];
        }

        if (iv != nsys - 1) {
            /* ==== eliminate the VB(iv+1) block, rows 1 -> 3 */
            for (int k = 0; k < 3; k++) {
                const double vtmp1 = xf->vb[ivp][0][k];
                const double vtmp2 = xf->vb[ivp][1][k];
                const double vtmp3 = vm[ivp][k][iv];
                row_sub3(vm[ivp][k], vm[iv][0], vm[iv][1], vm[iv][2], vtmp1, vtmp2, vtmp3, ivp, nsys);
                xf->vdel[ivp][0][k] = xf->vdel[ivp][0][k]
                                  - (vtmp1 * xf->vdel[iv][0][0] + vtmp2 * xf->vdel[iv][0][1] + vtmp3 * xf->vdel[iv][0][2]);
                xf->vdel[ivp][1][k] = xf->vdel[ivp][1][k]
                                  - (vtmp1 * xf->vdel[iv][1][0] + vtmp2 * xf->vdel[iv][1][1] + vtmp3 * xf->vdel[iv][1][2]);
            }

            if (iv == ivte1) {
                /* eliminate the VZ block */
                const int ivz = xf->isys[1][xf->iblte[1] + 1];
                for (int k = 0; k < 3; k++) {
                    const double vtmp1 = xf->vz[0][k];
                    const double vtmp2 = xf->vz[1][k];
                    row_sub2(vm[ivz][k], vm[iv][0], vm[iv][1], vtmp1, vtmp2, ivp, nsys);
                    xf->vdel[ivz][0][k] = xf->vdel[ivz][0][k] - (vtmp1 * xf->vdel[iv][0][0] + vtmp2 * xf->vdel[iv][0][1]);
                    xf->vdel[ivz][1][k] = xf->vdel[ivz][1][k] - (vtmp1 * xf->vdel[iv][1][0] + vtmp2 * xf->vdel[iv][1][1]);
                }
            }

            if (ivp != nsys - 1) {
                /* ==== eliminate the lower VM column */
                for (int kv = iv + 2; kv < nsys; kv++) {
                    const double vtmp1 = vm[kv][0][iv];
                    const double vtmp2 = vm[kv][1][iv];
                    const double vtmp3 = vm[kv][2][iv];

                    if (fabs(vtmp1) > vacc1) {
                        row_sub1(vm[kv][0], vm[iv][2], vtmp1, ivp, nsys);
                        xf->vdel[kv][0][0] = xf->vdel[kv][0][0] - vtmp1 * xf->vdel[iv][0][2];
                        xf->vdel[kv][1][0] = xf->vdel[kv][1][0] - vtmp1 * xf->vdel[iv][1][2];
                    }
                    if (fabs(vtmp2) > vacc2) {
                        row_sub1(vm[kv][1], vm[iv][2], vtmp2, ivp, nsys);
                        xf->vdel[kv][0][1] = xf->vdel[kv][0][1] - vtmp2 * xf->vdel[iv][0][2];
                        xf->vdel[kv][1][1] = xf->vdel[kv][1][1] - vtmp2 * xf->vdel[iv][1][2];
                    }
                    if (fabs(vtmp3) > vacc3) {
                        row_sub1(vm[kv][2], vm[iv][2], vtmp3, ivp, nsys);
                        xf->vdel[kv][0][2] = xf->vdel[kv][0][2] - vtmp3 * xf->vdel[iv][0][2];
                        xf->vdel[kv][1][2] = xf->vdel[kv][1][2] - vtmp3 * xf->vdel[iv][1][2];
                    }
                }
            }
        }
    }

    for (int iv = nsys - 1; iv >= 1; iv--) {
        /* eliminate the upper VM columns */
        double vtmp = xf->vdel[iv][0][2];
        for (int kv = iv - 1; kv >= 0; kv--) {
            xf->vdel[kv][0][0] = xf->vdel[kv][0][0] - vm[kv][0][iv] * vtmp;
            xf->vdel[kv][0][1] = xf->vdel[kv][0][1] - vm[kv][1][iv] * vtmp;
            xf->vdel[kv][0][2] = xf->vdel[kv][0][2] - vm[kv][2][iv] * vtmp;
        }

        vtmp = xf->vdel[iv][1][2];
        for (int kv = iv - 1; kv >= 0; kv--) {
            xf->vdel[kv][1][0] = xf->vdel[kv][1][0] - vm[kv][0][iv] * vtmp;
            xf->vdel[kv][1][1] = xf->vdel[kv][1][1] - vm[kv][1][iv] * vtmp;
            xf->vdel[kv][1][2] = xf->vdel[kv][1][2] - vm[kv][2][iv] * vtmp;
        }
    }
}
#else
void blsolv(xfoil_t *xf)
{
    /* Blocks, all 0-based: Fortran VA(k,l,iv) is va[iv][l-1][k-1] (row k, column l),
     * VM(k,jv,iv) is vm[iv][jv][k-1], VDEL(k,l,iv) is vdel[iv][l-1][k-1]. */

    const int ivte1 = xf->isys[0][xf->iblte[0]];

    const double vacc1 = xf->vaccel;
    const double vacc2 = xf->vaccel * 2.0 / (xf->s[xf->n - 1] - xf->s[0]);
    const double vacc3 = xf->vaccel * 2.0 / (xf->s[xf->n - 1] - xf->s[0]);

    for (int iv = 0; iv < xf->nsys; iv++) {
        const int ivp = iv + 1;

        /* ==== invert the VA(iv) block */

        /* normalize the first row */
        double pivot = 1.0 / xf->va[iv][0][0];
        xf->va[iv][1][0] = xf->va[iv][1][0] * pivot;
        for (int l = iv; l < xf->nsys; l++) xf->vm[iv][l][0] = xf->vm[iv][l][0] * pivot;
        xf->vdel[iv][0][0] = xf->vdel[iv][0][0] * pivot;
        xf->vdel[iv][1][0] = xf->vdel[iv][1][0] * pivot;

        /* eliminate the lower first column in the VA block */
        for (int k = 1; k < 3; k++) {
            const double vtmp = xf->va[iv][0][k];
            xf->va[iv][1][k] = xf->va[iv][1][k] - vtmp * xf->va[iv][1][0];
            for (int l = iv; l < xf->nsys; l++) xf->vm[iv][l][k] = xf->vm[iv][l][k] - vtmp * xf->vm[iv][l][0];
            xf->vdel[iv][0][k] = xf->vdel[iv][0][k] - vtmp * xf->vdel[iv][0][0];
            xf->vdel[iv][1][k] = xf->vdel[iv][1][k] - vtmp * xf->vdel[iv][1][0];
        }

        /* normalize the second row */
        pivot = 1.0 / xf->va[iv][1][1];
        for (int l = iv; l < xf->nsys; l++) xf->vm[iv][l][1] = xf->vm[iv][l][1] * pivot;
        xf->vdel[iv][0][1] = xf->vdel[iv][0][1] * pivot;
        xf->vdel[iv][1][1] = xf->vdel[iv][1][1] * pivot;

        /* eliminate the lower second column in the VA block */
        {
            const int k = 2;   /* third row */
            const double vtmp = xf->va[iv][1][k];
            for (int l = iv; l < xf->nsys; l++) xf->vm[iv][l][k] = xf->vm[iv][l][k] - vtmp * xf->vm[iv][l][1];
            xf->vdel[iv][0][k] = xf->vdel[iv][0][k] - vtmp * xf->vdel[iv][0][1];
            xf->vdel[iv][1][k] = xf->vdel[iv][1][k] - vtmp * xf->vdel[iv][1][1];
        }

        /* normalize the third row */
        pivot = 1.0 / xf->vm[iv][iv][2];
        for (int l = ivp; l < xf->nsys; l++) xf->vm[iv][l][2] = xf->vm[iv][l][2] * pivot;
        xf->vdel[iv][0][2] = xf->vdel[iv][0][2] * pivot;
        xf->vdel[iv][1][2] = xf->vdel[iv][1][2] * pivot;

        /* eliminate the upper third column in the VA block */
        {
            const double vtmp1 = xf->vm[iv][iv][0];
            const double vtmp2 = xf->vm[iv][iv][1];
            for (int l = ivp; l < xf->nsys; l++) {
                xf->vm[iv][l][0] = xf->vm[iv][l][0] - vtmp1 * xf->vm[iv][l][2];
                xf->vm[iv][l][1] = xf->vm[iv][l][1] - vtmp2 * xf->vm[iv][l][2];
            }
            xf->vdel[iv][0][0] = xf->vdel[iv][0][0] - vtmp1 * xf->vdel[iv][0][2];
            xf->vdel[iv][0][1] = xf->vdel[iv][0][1] - vtmp2 * xf->vdel[iv][0][2];
            xf->vdel[iv][1][0] = xf->vdel[iv][1][0] - vtmp1 * xf->vdel[iv][1][2];
            xf->vdel[iv][1][1] = xf->vdel[iv][1][1] - vtmp2 * xf->vdel[iv][1][2];
        }

        /* eliminate the upper second column in the VA block */
        {
            const double vtmp = xf->va[iv][1][0];
            for (int l = ivp; l < xf->nsys; l++) xf->vm[iv][l][0] = xf->vm[iv][l][0] - vtmp * xf->vm[iv][l][1];
            xf->vdel[iv][0][0] = xf->vdel[iv][0][0] - vtmp * xf->vdel[iv][0][1];
            xf->vdel[iv][1][0] = xf->vdel[iv][1][0] - vtmp * xf->vdel[iv][1][1];
        }

        if (iv != xf->nsys - 1) {
            /* ==== eliminate the VB(iv+1) block, rows 1 -> 3 */
            for (int k = 0; k < 3; k++) {
                const double vtmp1 = xf->vb[ivp][0][k];
                const double vtmp2 = xf->vb[ivp][1][k];
                const double vtmp3 = xf->vm[ivp][iv][k];
                for (int l = ivp; l < xf->nsys; l++)
                    xf->vm[ivp][l][k] = xf->vm[ivp][l][k] - (vtmp1 * xf->vm[iv][l][0] + vtmp2 * xf->vm[iv][l][1] + vtmp3 * xf->vm[iv][l][2]);
                xf->vdel[ivp][0][k] = xf->vdel[ivp][0][k]
                                  - (vtmp1 * xf->vdel[iv][0][0] + vtmp2 * xf->vdel[iv][0][1] + vtmp3 * xf->vdel[iv][0][2]);
                xf->vdel[ivp][1][k] = xf->vdel[ivp][1][k]
                                  - (vtmp1 * xf->vdel[iv][1][0] + vtmp2 * xf->vdel[iv][1][1] + vtmp3 * xf->vdel[iv][1][2]);
            }

            if (iv == ivte1) {
                /* eliminate the VZ block */
                const int ivz = xf->isys[1][xf->iblte[1] + 1];
                for (int k = 0; k < 3; k++) {
                    const double vtmp1 = xf->vz[0][k];
                    const double vtmp2 = xf->vz[1][k];
                    for (int l = ivp; l < xf->nsys; l++)
                        xf->vm[ivz][l][k] = xf->vm[ivz][l][k] - (vtmp1 * xf->vm[iv][l][0] + vtmp2 * xf->vm[iv][l][1]);
                    xf->vdel[ivz][0][k] = xf->vdel[ivz][0][k] - (vtmp1 * xf->vdel[iv][0][0] + vtmp2 * xf->vdel[iv][0][1]);
                    xf->vdel[ivz][1][k] = xf->vdel[ivz][1][k] - (vtmp1 * xf->vdel[iv][1][0] + vtmp2 * xf->vdel[iv][1][1]);
                }
            }

            if (ivp != xf->nsys - 1) {
                /* ==== eliminate the lower VM column */
                for (int kv = iv + 2; kv < xf->nsys; kv++) {
                    const double vtmp1 = xf->vm[kv][iv][0];
                    const double vtmp2 = xf->vm[kv][iv][1];
                    const double vtmp3 = xf->vm[kv][iv][2];

                    if (fabs(vtmp1) > vacc1) {
                        for (int l = ivp; l < xf->nsys; l++) xf->vm[kv][l][0] = xf->vm[kv][l][0] - vtmp1 * xf->vm[iv][l][2];
                        xf->vdel[kv][0][0] = xf->vdel[kv][0][0] - vtmp1 * xf->vdel[iv][0][2];
                        xf->vdel[kv][1][0] = xf->vdel[kv][1][0] - vtmp1 * xf->vdel[iv][1][2];
                    }
                    if (fabs(vtmp2) > vacc2) {
                        for (int l = ivp; l < xf->nsys; l++) xf->vm[kv][l][1] = xf->vm[kv][l][1] - vtmp2 * xf->vm[iv][l][2];
                        xf->vdel[kv][0][1] = xf->vdel[kv][0][1] - vtmp2 * xf->vdel[iv][0][2];
                        xf->vdel[kv][1][1] = xf->vdel[kv][1][1] - vtmp2 * xf->vdel[iv][1][2];
                    }
                    if (fabs(vtmp3) > vacc3) {
                        for (int l = ivp; l < xf->nsys; l++) xf->vm[kv][l][2] = xf->vm[kv][l][2] - vtmp3 * xf->vm[iv][l][2];
                        xf->vdel[kv][0][2] = xf->vdel[kv][0][2] - vtmp3 * xf->vdel[iv][0][2];
                        xf->vdel[kv][1][2] = xf->vdel[kv][1][2] - vtmp3 * xf->vdel[iv][1][2];
                    }
                }
            }
        }
    }

    for (int iv = xf->nsys - 1; iv >= 1; iv--) {
        /* eliminate the upper VM columns */
        double vtmp = xf->vdel[iv][0][2];
        for (int kv = iv - 1; kv >= 0; kv--) {
            xf->vdel[kv][0][0] = xf->vdel[kv][0][0] - xf->vm[kv][iv][0] * vtmp;
            xf->vdel[kv][0][1] = xf->vdel[kv][0][1] - xf->vm[kv][iv][1] * vtmp;
            xf->vdel[kv][0][2] = xf->vdel[kv][0][2] - xf->vm[kv][iv][2] * vtmp;
        }

        vtmp = xf->vdel[iv][1][2];
        for (int kv = iv - 1; kv >= 0; kv--) {
            xf->vdel[kv][1][0] = xf->vdel[kv][1][0] - xf->vm[kv][iv][0] * vtmp;
            xf->vdel[kv][1][1] = xf->vdel[kv][1][1] - xf->vm[kv][iv][1] * vtmp;
            xf->vdel[kv][1][2] = xf->vdel[kv][1][2] - xf->vm[kv][iv][2] * vtmp;
        }
    }
}
#endif
