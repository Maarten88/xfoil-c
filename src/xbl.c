/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of m_xbl.f90. Part 1: dslim, blpini, xifset.
 */
#include "xbl.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "xstate.h"
#include "xlayout.h"
#include "xblsys.h"
#include "xfortran.h"
#include "sxfoil.h"
#include "xpanel.h"
#include "xsolve.h"
#include "xspline.h"

#define P2(v) __builtin_powi((v), 2)
#define P3(v) __builtin_powi((v), 3)

void dslim(double *dstr, double thet, double uedg, double msq, double hklim)
{
    (void)uedg;
    double hk, hk_h, hk_m;
    const double h = *dstr / thet;
    hkin(h, msq, &hk, &hk_h, &hk_m);

    const double dh = fmax_f(0.0, hklim - hk) / hk_h;
    *dstr = *dstr + dh * thet;
}

void blpini(xfoil_t *xf)
{
    xf->sccon = 5.6;
    xf->gacon = 6.70;
    xf->gbcon = 0.75;
    xf->gccon = 18.0;
    xf->dlcon = 0.9;

    xf->ctrcon = 1.8;
    xf->ctrcex = 3.3;

    xf->duxcon = 1.0;

    xf->ctcon = 0.5 / (P2(xf->gacon) * xf->gbcon);

    xf->cffac = 1.0;
}

void xifset(xfoil_t *xf, int is)
{
    /* forced transition at (or beyond) the TE */
    if (xf->xstrip[is] >= 1.0) {
        xf->xiforc = xf->xssi[is][xf->iblte[is]];
        return;
    }

    const double chx = xf->xte - xf->xle;
    const double chy = xf->yte - xf->yle;
    const double chsq = P2(chx) + P2(chy);

    /* chord-based x/c, y/c */
    for (int i = 0; i < xf->n; i++) {
        xf->w1[i] = ((xf->x[i] - xf->xle) * chx + (xf->y[i] - xf->yle) * chy) / chsq;
        xf->w2[i] = ((xf->y[i] - xf->yle) * chx - (xf->x[i] - xf->xle) * chy) / chsq;
    }

    splind(xf->w1, xf->w3, xf->s, xf->n, -999.0, -999.0);
    splind(xf->w2, xf->w4, xf->s, xf->n, -999.0, -999.0);

    double str;
    if (is == 0) {
        /* approximate arc length of the forced transition point, then the actual one */
        str = xf->sle + (xf->s[0] - xf->sle) * xf->xstrip[is];
        sinvrt(&str, xf->xstrip[is], xf->w1, xf->w3, xf->s, xf->n);

        /* BL coordinate value */
        xf->xiforc = fmin_f((xf->sst - str), xf->xssi[is][xf->iblte[is]]);
    } else {
        /* same for the bottom side */
        str = xf->sle + (xf->s[xf->n - 1] - xf->sle) * xf->xstrip[is];
        sinvrt(&str, xf->xstrip[is], xf->w1, xf->w3, xf->s, xf->n);
        xf->xiforc = fmin_f((str - xf->sst), xf->xssi[is][xf->iblte[is]]);
    }

    if (xf->xiforc < 0.0) {
        if (xfoil_print) {
            printf("\n ***  Stagnation point is past trip on side%2d  ***\n", is + 1);
        }
        xf->xiforc = xf->xssi[is][xf->iblte[is]];
    }
}

/* ------------------------------------------------------------------------- */
/* Part 2: mrchue. Stations b are 0-based (Fortran ibl = b + 1); sides is 0/1.
 * Printed station and side numbers are 1-based, as in the Fortran. */

/* gfortran 3-argument MAX(a, b, c) is a left fold of the 2-argument one (probed) */
static double fmax3_f(double a, double b, double c)
{
    return fmax_f(fmax_f(a, b), c);
}

int mrchue(xfoil_t *xf)
{
    int result = 1;

    /* shape parameters for separation criteria */
    const double hlmax = 3.8;
    const double htmax = 2.5;

    /* Fortran locals live for the whole call: values carry over from side 1 to side 2 */
    double cte = 0.0, tte = 0.0, dte = 0.0, htarg = 0.0, hmax = 0.0, dmax = 0.0;

    for (int is = 0; is < 2; is++) {
        if (xfoil_print) printf(" %s%12d %s\n", "   side ", is + 1, " ...");   /* list-directed */

        xf->amcrit = xf->acrit[is];

        /* forced transition arc length position */
        xifset(xf, is);

        /* initialize the similarity station with Thwaites' formula */
        int b = 1;
        double xsi = xf->xssi[is][b];
        double uei = xf->uedg[is][b];
        xf->bule = 1.0;
        const double ucon = uei / pow(xsi, xf->bule);
        const double tsq = 0.45 / (ucon * (5.0 * xf->bule + 1.0) * xf->reybl) * pow(xsi, 1.0 - xf->bule);
        double thi = sqrt(tsq);
        double dsi = 2.2 * thi;
        double ami = 0.0;

        /* initial Ctau for the first turbulent station */
        double cti = 0.03;

        xf->tran = 0;
        xf->turb = 0;
        xf->itran[is] = xf->iblte[is];

        /* march downstream */
        for (b = 1; b < xf->nbl[is]; b++) {
            const int ibm = b - 1;

            xf->simi = (b == 1);
            xf->wake = (b > xf->iblte[is]);

            /* prescribed quantities */
            xsi = xf->xssi[is][b];
            uei = xf->uedg[is][b];

            double dswaki;
            if (xf->wake) {
                const int iw = b - xf->iblte[is];          /* Fortran IW = IBL - IBLTE (same value) */
                dswaki = xf->wgap[iw - 1];
            } else {
                dswaki = 0.0;
            }

            int direct = 1;

            /* Newton iteration loop for the current station */
            int converged = 0;
            for (int itbl = 1; itbl <= 25; itbl++) {
                /* assemble the 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi at the
                 * previous "1" and the current "2" station (the "1" coefficients are ignored) */
                blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(xf);

                /* check for transition and set the appropriate flags */
                if (!xf->simi && !xf->turb) {
                    result = trchek(xf);
                    if (ABORT_ON_NAN && !result) return result;
                    ami = xf->com2.ampl;

                    if (xf->tran) {
                        xf->itran[is] = b;
                        if (cti <= 0.0) {
                            cti = 0.03;
                            xf->com2.s = cti;
                        }
                    } else {
                        xf->itran[is] = b + 2;
                    }
                }

                if (b == xf->iblte[is] + 1) {
                    tte = xf->thet[0][xf->iblte[0]] + xf->thet[1][xf->iblte[1]];
                    dte = xf->dstr[0][xf->iblte[0]] + xf->dstr[1][xf->iblte[1]] + xf->ante;
                    cte = (xf->ctau[0][xf->iblte[0]] * xf->thet[0][xf->iblte[0]] + xf->ctau[1][xf->iblte[1]] * xf->thet[1][xf->iblte[1]]) / tte;
                    tesys(xf, cte, tte, dte);
                } else {
                    blsys(xf);
                }

                if (direct) {
                    /* try direct mode (dUe = 0 in the currently empty 4th line) */
                    xf->vs2[0][3] = 0.0;
                    xf->vs2[1][3] = 0.0;
                    xf->vs2[2][3] = 0.0;
                    xf->vs2[3][3] = 1.0;
                    xf->vsrez[3] = 0.0;

                    /* solve the Newton system for the current "2" station */
                    gauss(4, 4, &xf->vs2[0][0], xf->vsrez, 1);

                    /* max changes, underrelax if necessary */
                    dmax = fmax_f(fabs(xf->vsrez[1] / thi), fabs(xf->vsrez[2] / dsi));
                    if (b < xf->itran[is]) dmax = fmax_f(dmax, fabs(xf->vsrez[0] / 10.0));
                    if (b >= xf->itran[is]) dmax = fmax_f(dmax, fabs(xf->vsrez[0] / cti));

                    xf->rlx = 1.0;
                    if (dmax > 0.3) xf->rlx = 0.3 / dmax;

                    /* is direct mode not applicable? */
                    if (b != xf->iblte[is] + 1) {
                        /* resulting kinematic shape parameter Hk */
                        const double msq = uei * uei * xf->hstinv / (xf->gm1bl * (1.0 - 0.5 * uei * uei * xf->hstinv));
                        const double htest = (dsi + xf->rlx * xf->vsrez[2]) / (thi + xf->rlx * xf->vsrez[1]);
                        double hktest, dummy;
                        hkin(htest, msq, &hktest, &dummy, &dummy);

                        /* direct or inverse problem, based on Hk */
                        if (b < xf->itran[is]) hmax = hlmax;
                        if (b >= xf->itran[is]) hmax = htmax;
                        direct = hktest < hmax;
                    }

                    if (direct) {
                        /* update as usual */
                        if (b >= xf->itran[is]) cti = cti + xf->rlx * xf->vsrez[0];
                        thi = thi + xf->rlx * xf->vsrez[1];
                        dsi = dsi + xf->rlx * xf->vsrez[2];
                    } else {
                        /* prescribed Hk for the inverse calculation at the current station */
                        if (b < xf->itran[is]) {
                            /* laminar: relatively slow increase in Hk downstream */
                            htarg = xf->com1.hk + 0.03 * (xf->com2.x - xf->com1.x) / xf->com1.t;
                        } else if (b == xf->itran[is]) {
                            /* transition interval: weighted laminar and turbulent case */
                            htarg = xf->com1.hk + (0.03 * (xf->xt - xf->com1.x) - 0.15 * (xf->com2.x - xf->xt)) / xf->com1.t;
                        } else if (xf->wake) {
                            /* turbulent wake: asymptotic wake behavior, approximate backward Euler */
                            const double cst = 0.03 * (xf->com2.x - xf->com1.x) / xf->com1.t;
                            xf->com2.hk = xf->com1.hk;
                            xf->com2.hk = xf->com2.hk - (xf->com2.hk + cst * P3(xf->com2.hk - 1.0) - xf->com1.hk)
                                                / (1.0 + 3.0 * cst * P2(xf->com2.hk - 1.0));
                            xf->com2.hk = xf->com2.hk - (xf->com2.hk + cst * P3(xf->com2.hk - 1.0) - xf->com1.hk)
                                                / (1.0 + 3.0 * cst * P2(xf->com2.hk - 1.0));
                            xf->com2.hk = xf->com2.hk - (xf->com2.hk + cst * P3(xf->com2.hk - 1.0) - xf->com1.hk)
                                                / (1.0 + 3.0 * cst * P2(xf->com2.hk - 1.0));
                            htarg = xf->com2.hk;
                        } else {
                            /* turbulent: relatively fast decrease in Hk downstream */
                            htarg = xf->com1.hk - 0.15 * (xf->com2.x - xf->com1.x) / xf->com1.t;
                        }

                        /* limit the specified Hk to something reasonable */
                        if (xf->wake)
                            htarg = fmax_f(htarg, 1.01);
                        else
                            htarg = fmax_f(htarg, hmax);

                        if (xfoil_print) {
                            printf(" MRCHUE: Inverse mode at%4d     Hk =%8.3f\n", b + 1, htarg);
                        }

                        /* try again with the prescribed Hk */
                        continue;
                    }
                } else {
                    /* inverse mode (force Hk to the prescribed value HTARG) */
                    xf->vs2[0][3] = 0.0;
                    xf->vs2[1][3] = xf->com2.hk_t;
                    xf->vs2[2][3] = xf->com2.hk_d;
                    xf->vs2[3][3] = xf->com2.hk_u;
                    xf->vsrez[3] = htarg - xf->com2.hk;

                    gauss(4, 4, &xf->vs2[0][0], xf->vsrez, 1);

                    /* added Ue clamp, MD 3 Apr 03 */
                    dmax = fmax3_f(fabs(xf->vsrez[1] / thi), fabs(xf->vsrez[2] / dsi), fabs(xf->vsrez[3] / uei));
                    if (b >= xf->itran[is]) dmax = fmax_f(dmax, fabs(xf->vsrez[0] / cti));

                    xf->rlx = 1.0;
                    if (dmax > 0.3) xf->rlx = 0.3 / dmax;

                    /* update variables */
                    if (b >= xf->itran[is]) cti = cti + xf->rlx * xf->vsrez[0];
                    thi = thi + xf->rlx * xf->vsrez[1];
                    dsi = dsi + xf->rlx * xf->vsrez[2];
                    uei = uei + xf->rlx * xf->vsrez[3];
                }

                /* eliminate absurd transients */
                if (b >= xf->itran[is]) {
                    /* Fortran min(cti, 0.30): gfortran gives no NaN guarantee for MIN/MAX; this
                     * site compiles to `minsd cti, 0.30` = (0.30 < cti) ? 0.30 : cti, which keeps a
                     * NaN cti (verified with objdump; diverged-march states depend on it) */
                    cti = (0.30 < cti) ? 0.30 : cti;
                    cti = fmax_f(cti, 0.0000001);
                }

                const double hklim = (b <= xf->iblte[is]) ? 1.02 : 1.00005;
                const double msq = uei * uei * xf->hstinv / (xf->gm1bl * (1.0 - 0.5 * uei * uei * xf->hstinv));
                double dsw = dsi - dswaki;
                dslim(&dsw, thi, uei, msq, hklim);
                dsi = dsw + dswaki;

                if (dmax <= 1.0e-5) {
                    converged = 1;
                    break;
                }
            }

            if (!converged) {
                if (xfoil_print) {
                    printf(" MRCHUE: Convergence failed at%4d  side%2d    Res =%12.4E\n", b + 1, is + 1, dmax);
                }

                /* the current unconverged solution might still be reasonable... */
                if (dmax > 0.1) {
                    /* ...but it is garbage: extrapolate values instead */
                    if (b > 2) {
                        if (b <= xf->iblte[is]) {
                            thi = xf->thet[is][ibm] * pow(xf->xssi[is][b] / xf->xssi[is][ibm], 0.5);
                            dsi = xf->dstr[is][ibm] * pow(xf->xssi[is][b] / xf->xssi[is][ibm], 0.5);
                        } else if (b == xf->iblte[is] + 1) {
                            cti = cte;
                            thi = tte;
                            dsi = dte;
                        } else {
                            thi = xf->thet[is][ibm];
                            const double ratlen = (xf->xssi[is][b] - xf->xssi[is][ibm]) / (10.0 * xf->dstr[is][ibm]);
                            dsi = (xf->dstr[is][ibm] + thi * ratlen) / (1.0 + ratlen);
                        }
                        if (b == xf->itran[is]) cti = 0.05;
                        if (b > xf->itran[is]) cti = xf->ctau[is][ibm];

                        uei = xf->uedg[is][b];
                        if (b > 1 && b < xf->nbl[is] - 1) uei = 0.5 * (xf->uedg[is][b - 1] + xf->uedg[is][b + 1]);
                    }
                }

                blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(xf);

                /* check for transition and set the appropriate flags */
                if (!xf->simi && !xf->turb) {
                    result = trchek(xf);
                    if (ABORT_ON_NAN && !result) return result;
                    ami = xf->com2.ampl;
                    if (xf->tran) xf->itran[is] = b;
                    if (!xf->tran) xf->itran[is] = b + 2;
                }

                /* all other extrapolated values for the current station */
                if (b < xf->itran[is]) blvar(xf, 1);
                if (b >= xf->itran[is]) blvar(xf, 2);
                if (xf->wake) blvar(xf, 3);

                if (b < xf->itran[is]) blmid(xf, 1);
                if (b >= xf->itran[is]) blmid(xf, 2);
                if (xf->wake) blmid(xf, 3);
            }

            /* store primary variables (after the Newton iterations) */
            if (b < xf->itran[is]) xf->ctau[is][b] = ami;
            if (b >= xf->itran[is]) xf->ctau[is][b] = cti;
            xf->thet[is][b] = thi;
            xf->dstr[is][b] = dsi;
            xf->uedg[is][b] = uei;
            xf->mass[is][b] = dsi * uei;
            xf->tau[is][b] = 0.5 * xf->com2.r * xf->com2.u * xf->com2.u * xf->com2.cf;
            xf->dis[is][b] = xf->com2.r * xf->com2.u * xf->com2.u * xf->com2.u * xf->com2.di * xf->com2.hs * 0.5;
            xf->ctq[is][b] = xf->com2.cq;
            xf->delt[is][b] = xf->com2.de;
            xf->tstr[is][b] = xf->com2.hs * xf->com2.t;

            /* "1" variables := "2" variables for the next streamwise station */
            blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
            blkin(xf);
            xf->com1 = xf->com2;

            /* turbulent intervals follow the transition interval or the TE */
            if (xf->tran || b == xf->iblte[is]) {
                xf->turb = 1;
                /* save the transition location */
                xf->tforce[is] = xf->trforc;
            }

            xf->tran = 0;

            if (b == xf->iblte[is]) {
                thi = xf->thet[0][xf->iblte[0]] + xf->thet[1][xf->iblte[1]];
                dsi = xf->dstr[0][xf->iblte[0]] + xf->dstr[1][xf->iblte[1]] + xf->ante;
            }
        }
    }

    return result;
}

/* ------------------------------------------------------------------------- */
/* Part 3: mrchdu. Its min/max sites compile (objdump of gfortran's mrchdu) to the
 * standard orders fmin_f(a,b) = (a<b)?a:b and fmax_f(a,b) = (a>b)?a:b, unlike the
 * cti clamp in mrchue. */

int mrchdu(xfoil_t *xf)
{
    static const double deps = 5.0e-6;

    int result = 1;
    double vtmp[5][4], vztmp[4];

    /* constant controlling how far Hk may deviate from the specified value */
    const double senswt = 1000.0;
    double sens = 0.0, sennew = 0.0;

    /* Fortran locals live for the whole call: values carry over from side 1 to side 2 */
    double ami = 0.0, cti = 0.0, cte = 0.0, tte = 0.0, dte = 0.0, dmax = 0.0;
    double ueref = 0.0, hkref = 0.0;

    for (int is = 0; is < 2; is++) {
        xf->amcrit = xf->acrit[is];

        /* forced transition arc length position */
        xifset(xf, is);

        /* leading edge pressure gradient parameter x/u du/dx */
        int b = 1;
        double xsi = xf->xssi[is][b];
        double uei = xf->uedg[is][b];
        xf->bule = 1.0;

        /* old transition station */
        const int itrold = xf->itran[is];

        xf->tran = 0;
        xf->turb = 0;
        xf->itran[is] = xf->iblte[is];

        /* march downstream */
        for (b = 1; b < xf->nbl[is]; b++) {
            const int ibm = b - 1;

            xf->simi = (b == 1);
            xf->wake = (b > xf->iblte[is]);

            /* initialize the current station to the existing variables */
            xsi = xf->xssi[is][b];
            uei = xf->uedg[is][b];
            double thi = xf->thet[is][b];
            double dsi = xf->dstr[is][b];

            /* fixed BUG, MD 7 June 99 */
            if (b < itrold) {
                ami = xf->ctau[is][b];
                cti = 0.03;
            } else {
                cti = xf->ctau[is][b];
                if (cti <= 0.0) cti = 0.03;
            }

            double dswaki;
            if (xf->wake) {
                const int iw = b - xf->iblte[is];          /* Fortran IW = IBL - IBLTE (same value) */
                dswaki = xf->wgap[iw - 1];
            } else {
                dswaki = 0.0;
            }

            if (b <= xf->iblte[is]) dsi = fmax_f(dsi - dswaki, 1.02000 * thi) + dswaki;
            if (b > xf->iblte[is]) dsi = fmax_f(dsi - dswaki, 1.00005 * thi) + dswaki;

            /* Newton iteration loop for the current station */
            int converged = 0;
            for (int itbl = 1; itbl <= 25; itbl++) {
                /* assemble the 10x3 linearized system at the "1" and "2" stations */
                blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(xf);

                /* check for transition and set the appropriate flags */
                if (!xf->simi && !xf->turb) {
                    result = trchek(xf);
                    if (ABORT_ON_NAN && !result) return result;
                    ami = xf->com2.ampl;
                    if (xf->tran) xf->itran[is] = b;
                    if (!xf->tran) xf->itran[is] = b + 2;
                }

                if (b == xf->iblte[is] + 1) {
                    tte = xf->thet[0][xf->iblte[0]] + xf->thet[1][xf->iblte[1]];
                    dte = xf->dstr[0][xf->iblte[0]] + xf->dstr[1][xf->iblte[1]] + xf->ante;
                    cte = (xf->ctau[0][xf->iblte[0]] * xf->thet[0][xf->iblte[0]] + xf->ctau[1][xf->iblte[1]] * xf->thet[1][xf->iblte[1]]) / tte;
                    tesys(xf, cte, tte, dte);
                } else {
                    blsys(xf);
                }

                /* first iteration: "baseline" Ue and Hk for the Ue(Hk) relation */
                if (itbl == 1) {
                    ueref = xf->com2.u;
                    hkref = xf->com2.hk;

                    /* point was turbulent and is now laminar: extrapolate the baseline Hk */
                    if (b < xf->itran[is] && b >= itrold) {
                        const double uem = xf->uedg[is][b - 1];
                        const double dsm = xf->dstr[is][b - 1];
                        const double thm = xf->thet[is][b - 1];
                        const double msq = uem * uem * xf->hstinv / (xf->gm1bl * (1.0 - 0.5 * uem * uem * xf->hstinv));
                        double dummy;
                        hkin(dsm / thm, msq, &hkref, &dummy, &dummy);
                    }

                    /* point was laminar: reinitialize or extrapolate Ctau if it is now turbulent */
                    if (b < itrold) {
                        if (xf->tran) xf->ctau[is][b] = 0.03;
                        if (xf->turb) xf->ctau[is][b] = xf->ctau[is][b - 1];
                        if (xf->tran || xf->turb) {
                            cti = xf->ctau[is][b];
                            xf->com2.s = cti;
                        }
                    }
                }

                if (xf->simi || b == xf->iblte[is] + 1) {
                    /* similarity station or first wake point: prescribe Ue */
                    xf->vs2[0][3] = 0.0;
                    xf->vs2[1][3] = 0.0;
                    xf->vs2[2][3] = 0.0;
                    xf->vs2[3][3] = xf->com2.u_uei;
                    xf->vsrez[3] = ueref - xf->com2.u;
                } else {
                    /* Ue-Hk characteristic slope */
                    for (int k = 0; k < 4; k++) {
                        vztmp[k] = xf->vsrez[k];
                        for (int l = 0; l < 5; l++) vtmp[l][k] = xf->vs2[l][k];
                    }

                    /* unit dHk */
                    vtmp[0][3] = 0.0;
                    vtmp[1][3] = xf->com2.hk_t;
                    vtmp[2][3] = xf->com2.hk_d;
                    vtmp[3][3] = xf->com2.hk_u * xf->com2.u_uei;
                    vztmp[3] = 1.0;

                    /* dUe response */
                    gauss(4, 4, &vtmp[0][0], vztmp, 1);

                    /* SENSWT * (normalized dUe/dHk) */
                    sennew = senswt * vztmp[3] * hkref / ueref;
                    if (itbl <= 5)
                        sens = sennew;
                    else if (itbl <= 15)
                        sens = 0.5 * (sens + sennew);

                    /* prescribed Ue-Hk combination */
                    xf->vs2[0][3] = 0.0;
                    xf->vs2[1][3] = xf->com2.hk_t * hkref;
                    xf->vs2[2][3] = xf->com2.hk_d * hkref;
                    xf->vs2[3][3] = (xf->com2.hk_u * hkref + sens / ueref) * xf->com2.u_uei;
                    xf->vsrez[3] = -(P2(hkref) * (xf->com2.hk / hkref - 1.0)) - sens * (xf->com2.u / ueref - 1.0);
                }

                /* solve the Newton system for the current "2" station */
                gauss(4, 4, &xf->vs2[0][0], xf->vsrez, 1);

                /* max changes, underrelax if necessary (added Ue clamp, MD 3 Apr 03) */
                dmax = fmax3_f(fabs(xf->vsrez[1] / thi), fabs(xf->vsrez[2] / dsi), fabs(xf->vsrez[3] / uei));
                if (b >= xf->itran[is]) dmax = fmax_f(dmax, fabs(xf->vsrez[0] / (10.0 * cti)));

                xf->rlx = 1.0;
                if (dmax > 0.3) xf->rlx = 0.3 / dmax;

                /* update as usual */
                if (b < xf->itran[is]) ami = ami + xf->rlx * xf->vsrez[0];
                if (b >= xf->itran[is]) cti = cti + xf->rlx * xf->vsrez[0];
                thi = thi + xf->rlx * xf->vsrez[1];
                dsi = dsi + xf->rlx * xf->vsrez[2];
                uei = uei + xf->rlx * xf->vsrez[3];

                /* eliminate absurd transients */
                if (b >= xf->itran[is]) {
                    cti = fmin_f(cti, 0.30);
                    cti = fmax_f(cti, 0.0000001);
                }

                const double hklim = (b <= xf->iblte[is]) ? 1.02 : 1.00005;
                const double msq = uei * uei * xf->hstinv / (xf->gm1bl * (1.0 - 0.5 * uei * uei * xf->hstinv));
                double dsw = dsi - dswaki;
                dslim(&dsw, thi, uei, msq, hklim);
                dsi = dsw + dswaki;

                if (dmax <= deps) {
                    converged = 1;
                    break;
                }
            }

            if (!converged) {
                if (xfoil_print) {
                    printf(" MRCHDU: Convergence failed at%4d  side%2d    Res =%12.4E\n", b + 1, is + 1, dmax);
                }

                /* the current unconverged solution might still be reasonable... */
                if (dmax > 0.1) {
                    /* ...but it is garbage: extrapolate values instead */
                    if (b > 2) {
                        if (b <= xf->iblte[is]) {
                            thi = xf->thet[is][ibm] * pow(xf->xssi[is][b] / xf->xssi[is][ibm], 0.5);
                            dsi = xf->dstr[is][ibm] * pow(xf->xssi[is][b] / xf->xssi[is][ibm], 0.5);
                            uei = xf->uedg[is][ibm];
                        } else if (b == xf->iblte[is] + 1) {
                            cti = cte;
                            thi = tte;
                            dsi = dte;
                            uei = xf->uedg[is][ibm];
                        } else {
                            thi = xf->thet[is][ibm];
                            const double ratlen = (xf->xssi[is][b] - xf->xssi[is][ibm]) / (10.0 * xf->dstr[is][ibm]);
                            dsi = (xf->dstr[is][ibm] + thi * ratlen) / (1.0 + ratlen);
                            uei = xf->uedg[is][ibm];
                        }
                        if (b == xf->itran[is]) cti = 0.05;
                        if (b > xf->itran[is]) cti = xf->ctau[is][ibm];
                    }
                }

                blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(xf);

                /* check for transition and set the appropriate flags */
                if (!xf->simi && !xf->turb) {
                    result = trchek(xf);
                    if (ABORT_ON_NAN && !result) return result;
                    ami = xf->com2.ampl;
                    if (xf->tran) xf->itran[is] = b;
                    if (!xf->tran) xf->itran[is] = b + 2;
                }

                /* all other extrapolated values for the current station */
                if (b < xf->itran[is]) blvar(xf, 1);
                if (b >= xf->itran[is]) blvar(xf, 2);
                if (xf->wake) blvar(xf, 3);

                if (b < xf->itran[is]) blmid(xf, 1);
                if (b >= xf->itran[is]) blmid(xf, 2);
                if (xf->wake) blmid(xf, 3);
            }
            /* pick up here after the Newton iterations */
            sens = sennew;

            /* store primary variables */
            if (b < xf->itran[is]) xf->ctau[is][b] = ami;
            if (b >= xf->itran[is]) xf->ctau[is][b] = cti;
            xf->thet[is][b] = thi;
            xf->dstr[is][b] = dsi;
            xf->uedg[is][b] = uei;
            xf->mass[is][b] = dsi * uei;
            xf->tau[is][b] = 0.5 * xf->com2.r * xf->com2.u * xf->com2.u * xf->com2.cf;
            xf->dis[is][b] = xf->com2.r * xf->com2.u * xf->com2.u * xf->com2.u * xf->com2.di * xf->com2.hs * 0.5;
            xf->ctq[is][b] = xf->com2.cq;
            xf->delt[is][b] = xf->com2.de;
            xf->tstr[is][b] = xf->com2.hs * xf->com2.t;

            /* "1" variables := "2" variables for the next streamwise station */
            blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
            blkin(xf);
            xf->com1 = xf->com2;

            /* turbulent intervals follow the transition interval or the TE */
            if (xf->tran || b == xf->iblte[is]) {
                xf->turb = 1;
                /* save the transition location */
                xf->tforce[is] = xf->trforc;
            }

            xf->tran = 0;
        }
    }

    return result;
}

/* ------------------------------------------------------------------------- */
/* Part 4: setbl. Block storage as in xsolve.c:blsolv (va[iv][l][k], vm[iv][jv][k], ...,
 * all 0-based) and xblsys.c (vs1[l][k]). */

int setbl(xfoil_t *xf)
{
    /* work arrays in the context (xscratch.h) */
    double *d1_m = xf->scr.d1_m, *d2_m = xf->scr.d2_m, *u1_m = xf->scr.u1_m, *u2_m = xf->scr.u2_m;
    double *ule1_m = xf->scr.ule1_m, *ule2_m = xf->scr.ule2_m, *ute1_m = xf->scr.ute1_m, *ute2_m = xf->scr.ute2_m;
    double (*usav)[IVX] = xf->scr.usav;

    int result = 1;

    /* CL used to define the Mach and Reynolds numbers */
    const double clmr = xf->lalfa ? xf->cl : xf->clspec;

    /* current MINF(CL) */
    double ma_clmr, re_clmr;
    mrcl(xf, clmr, &ma_clmr, &re_clmr);
    const double msq_clmr = 2.0 * xf->minf * ma_clmr;

    /* compressibility parameter TKLAM and its derivative */
    comset(xf);

    /* gas constant (= Cp/Cv) */
    xf->gm1bl = xf->gamm1;

    /* parameters for the compressibility correction */
    xf->qinfbl = xf->qinf;
    xf->tkbl = xf->tklam;
    xf->tkbl_ms = xf->tkl_msq;

    /* stagnation density and 1/enthalpy */
    xf->rstbl = pow(1.0 + 0.5 * xf->gm1bl * P2(xf->minf), 1.0 / xf->gm1bl);
    xf->rstbl_ms = 0.5 * xf->rstbl / (1.0 + 0.5 * xf->gm1bl * P2(xf->minf));

    xf->hstinv = xf->gm1bl * P2(xf->minf / xf->qinfbl) / (1.0 + 0.5 * xf->gm1bl * P2(xf->minf));
    xf->hstinv_ms = xf->gm1bl * P2(1.0 / xf->qinfbl) / (1.0 + 0.5 * xf->gm1bl * P2(xf->minf))
                - 0.5 * xf->gm1bl * xf->hstinv / (1.0 + 0.5 * xf->gm1bl * P2(xf->minf));

    /* Reynolds number based on freestream density, velocity, viscosity */
    const double herat = 1.0 - 0.5 * P2(xf->qinfbl) * xf->hstinv;
    const double herat_ms = -(0.5 * P2(xf->qinfbl) * xf->hstinv_ms);

    xf->reybl = xf->reinf * sqrt(P3(herat)) * (1.0 + xf->hvrat) / (herat + xf->hvrat);
    xf->reybl_re = sqrt(P3(herat)) * (1.0 + xf->hvrat) / (herat + xf->hvrat);
    xf->reybl_ms = xf->reybl * (1.5 / herat - 1.0 / (herat + xf->hvrat)) * herat_ms;

    xf->idampv = xf->idamp;

    /* TE thickness */

    if (!xf->lblini) {
        /* initialize the BL by marching with Ue (fudge at separation) */
        if (xfoil_print) printf("\n Initializing BL ...\n");
        result = mrchue(xf);
        if (result)
            xf->lblini = 1;
        else if (ABORT_ON_NAN)
            return result;
    }

    if (xfoil_print) printf("\n");

    /* march the BL with the current Ue and Ds to establish transition */
    result = mrchdu(xf);
    if (ABORT_ON_NAN && !result) return result;

    for (int is = 0; is < 2; is++)
        for (int b = 1; b < xf->nbl[is]; b++) usav[is][b] = xf->uedg[is][b];

    ueset(xf);

    for (int is = 0; is < 2; is++) {
        for (int b = 1; b < xf->nbl[is]; b++) {
            const double temp = usav[is][b];
            usav[is][b] = xf->uedg[is][b];
            xf->uedg[is][b] = temp;
        }
    }

    const int ile1 = xf->ipan[0][1];
    const int ile2 = xf->ipan[1][1];
    const int ite1 = xf->ipan[0][xf->iblte[0]];
    const int ite2 = xf->ipan[1][xf->iblte[1]];

    const int jvte1 = xf->isys[0][xf->iblte[0]];
    const int jvte2 = xf->isys[1][xf->iblte[1]];

    const double dule1 = xf->uedg[0][1] - usav[0][1];
    const double dule2 = xf->uedg[1][1] - usav[1][1];

    /* LE and TE Ue sensitivities wrt all m values */
    for (int js = 0; js < 2; js++) {
        for (int jbl = 1; jbl < xf->nbl[js]; jbl++) {
            const int j = xf->ipan[js][jbl];
            const int jv = xf->isys[js][jbl];
            ule1_m[jv] = -(xf->vti[0][1] * xf->vti[js][jbl] * xf->dij[j][ile1]);
            ule2_m[jv] = -(xf->vti[1][1] * xf->vti[js][jbl] * xf->dij[j][ile2]);
            ute1_m[jv] = -(xf->vti[0][xf->iblte[0]] * xf->vti[js][jbl] * xf->dij[j][ite1]);
            ute2_m[jv] = -(xf->vti[1][xf->iblte[1]] * xf->vti[js][jbl] * xf->dij[j][ite2]);
        }
    }

    const double ule1_a = xf->uinv_a[0][1];
    const double ule2_a = xf->uinv_a[1][1];


    /* Fortran locals live for the whole call */
    double ami = 0.0, cti = 0.0;
    double tte = 0.0, dte = 0.0, cte = 0.0;
    double tte_tte1 = 0.0, tte_tte2 = 0.0, cte_cte1 = 0.0, cte_cte2 = 0.0, cte_tte1 = 0.0, cte_tte2 = 0.0;

    /* **** go over each boundary layer / wake */
    for (int is = 0; is < 2; is++) {
        /* there is no station "1" at similarity, so zero everything out */
        for (int js = 0; js < 2; js++) {
            for (int jbl = 1; jbl < xf->nbl[js]; jbl++) {
                const int jv = xf->isys[js][jbl];
                u1_m[jv] = 0.0;
                d1_m[jv] = 0.0;
            }
        }
        double u1_a = 0.0;
        double d1_a = 0.0;

        double due1 = 0.0;
        double dds1 = 0.0;

        /* similarity station pressure gradient parameter x/u du/dx */
        xf->bule = 1.0;

        xf->amcrit = xf->acrit[is];

        /* forced transition arc length position */
        xifset(xf, is);

        xf->tran = 0;
        xf->turb = 0;

        /* **** sweep downstream setting up the BL equation linearizations */
        for (int b = 1; b < xf->nbl[is]; b++) {
            const int iv = xf->isys[is][b];

            xf->simi = (b == 1);
            xf->wake = (b > xf->iblte[is]);
            xf->tran = (b == xf->itran[is]);
            xf->turb = (b > xf->itran[is]);

            const int i = xf->ipan[is][b];

            /* primary variables for the current station */
            const double xsi = xf->xssi[is][b];
            if (b < xf->itran[is]) ami = xf->ctau[is][b];
            if (b >= xf->itran[is]) cti = xf->ctau[is][b];
            const double uei = xf->uedg[is][b];
            const double thi = xf->thet[is][b];
            const double mdi = xf->mass[is][b];

            const double dsi = mdi / uei;

            double dswaki;
            if (xf->wake) {
                const int iw = b - xf->iblte[is];          /* Fortran IW = IBL - IBLTE (same value) */
                dswaki = xf->wgap[iw - 1];
            } else {
                dswaki = 0.0;
            }

            /* derivatives of DSI (= D2) */
            const double d2_m2 = 1.0 / uei;
            const double d2_u2 = -(dsi / uei);

            for (int js = 0; js < 2; js++) {
                for (int jbl = 1; jbl < xf->nbl[js]; jbl++) {
                    const int j = xf->ipan[js][jbl];
                    const int jv = xf->isys[js][jbl];
                    u2_m[jv] = -(xf->vti[is][b] * xf->vti[js][jbl] * xf->dij[j][i]);
                    d2_m[jv] = d2_u2 * u2_m[jv];
                }
            }
            d2_m[iv] = d2_m[iv] + d2_m2;

            const double u2_a = xf->uinv_a[is][b];
            const double d2_a = d2_u2 * u2_a;

            /* "forced" changes due to the mismatch between UEDG and USAV = UINV + dij*MASS */
            const double due2 = xf->uedg[is][b] - usav[is][b];
            const double dds2 = d2_u2 * due2;

            blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
            blkin(xf);

            /* check for transition and set TRAN, XT, etc. if found */
            if (xf->tran) {
                result = trchek(xf);
                if (ABORT_ON_NAN && !result) return result;
                ami = xf->com2.ampl;
            }
            if (b == xf->itran[is] && !xf->tran && xfoil_print) {
                printf(" SETBL: Xtr???  n1 n2:  %.17g %.17g\n", xf->com1.ampl, xf->com2.ampl);
            }

            /* assemble the 10x4 linearized system at the "1" and "2" stations */
            if (b == xf->iblte[is] + 1) {
                /* start of the wake: add the TE base thickness to Dstar */
                tte = xf->thet[0][xf->iblte[0]] + xf->thet[1][xf->iblte[1]];
                dte = xf->dstr[0][xf->iblte[0]] + xf->dstr[1][xf->iblte[1]] + xf->ante;
                cte = (xf->ctau[0][xf->iblte[0]] * xf->thet[0][xf->iblte[0]] + xf->ctau[1][xf->iblte[1]] * xf->thet[1][xf->iblte[1]]) / tte;
                tesys(xf, cte, tte, dte);

                tte_tte1 = 1.0;
                tte_tte2 = 1.0;
                const double dte_mte1 = 1.0 / xf->uedg[0][xf->iblte[0]];
                const double dte_ute1 = -(xf->dstr[0][xf->iblte[0]] / xf->uedg[0][xf->iblte[0]]);
                const double dte_mte2 = 1.0 / xf->uedg[1][xf->iblte[1]];
                const double dte_ute2 = -(xf->dstr[1][xf->iblte[1]] / xf->uedg[1][xf->iblte[1]]);
                cte_cte1 = xf->thet[0][xf->iblte[0]] / tte;
                cte_cte2 = xf->thet[1][xf->iblte[1]] / tte;
                cte_tte1 = (xf->ctau[0][xf->iblte[0]] - cte) / tte;
                cte_tte2 = (xf->ctau[1][xf->iblte[1]] - cte) / tte;

                /* re-define D1 sensitivities wrt m (D1 depends on both TE Ds values) */
                for (int js = 0; js < 2; js++) {
                    for (int jbl = 1; jbl < xf->nbl[js]; jbl++) {
                        const int jv = xf->isys[js][jbl];
                        d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
                    }
                }
                d1_m[jvte1] = d1_m[jvte1] + dte_mte1;
                d1_m[jvte2] = d1_m[jvte2] + dte_mte2;

                /* "forced" changes from the UEDG --- USAV mismatch */
                due1 = 0.0;
                dds1 = dte_ute1 * (xf->uedg[0][xf->iblte[0]] - usav[0][xf->iblte[0]])
                       + dte_ute2 * (xf->uedg[1][xf->iblte[1]] - usav[1][xf->iblte[1]]);
            } else {
                blsys(xf);
            }

            /* wall shear and equilibrium max shear coefficient for plotting output */
            xf->tau[is][b] = 0.5 * xf->com2.r * xf->com2.u * xf->com2.u * xf->com2.cf;
            xf->dis[is][b] = xf->com2.r * xf->com2.u * xf->com2.u * xf->com2.u * xf->com2.di * xf->com2.hs * 0.5;
            xf->ctq[is][b] = xf->com2.cq;
            xf->delt[is][b] = xf->com2.de;

            /* XI sensitivities wrt LE Ue changes */
            double xi_ule1, xi_ule2;
            if (is == 0) {
                xi_ule1 = xf->sst_go;
                xi_ule2 = -xf->sst_gp;
            } else {
                xi_ule1 = -xf->sst_go;
                xi_ule2 = xf->sst_gp;
            }

            /* stuff the BL system coefficients into the main Jacobian matrix */
            for (int k = 0; k < 3; k++) {
                /* coefficients as locals: the stores into VM could otherwise alias them (and nsys),
                 * forcing reloads in every iteration (same arithmetic) */
                const double c_d1 = xf->vs1[2][k], c_u1 = xf->vs1[3][k];
                const double c_d2 = xf->vs2[2][k], c_u2 = xf->vs2[3][k];
                const double c_x = xf->vs1[4][k] + xf->vs2[4][k] + xf->vsx[k];
                const int nsys = xf->nsys;
                for (int jv = 0; jv < nsys; jv++) {
                    XF_VM(xf, k, jv, iv) = c_d1 * d1_m[jv] + c_u1 * u1_m[jv]
                                    + c_d2 * d2_m[jv] + c_u2 * u2_m[jv]
                                    + c_x * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
                }

                xf->vb[iv][0][k] = xf->vs1[0][k];
                xf->vb[iv][1][k] = xf->vs1[1][k];

                xf->va[iv][0][k] = xf->vs2[0][k];
                xf->va[iv][1][k] = xf->vs2[1][k];

                if (xf->lalfa) {
                    xf->vdel[iv][1][k] = xf->vsr[k] * re_clmr + xf->vsm[k] * msq_clmr;
                } else {
                    xf->vdel[iv][1][k] = (xf->vs1[3][k] * u1_a + xf->vs1[2][k] * d1_a)
                                     + (xf->vs2[3][k] * u2_a + xf->vs2[2][k] * d2_a)
                                     + (xf->vs1[4][k] + xf->vs2[4][k] + xf->vsx[k]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
                }

                xf->vdel[iv][0][k] = xf->vsrez[k] + (xf->vs1[3][k] * due1 + xf->vs1[2][k] * dds1)
                                 + (xf->vs2[3][k] * due2 + xf->vs2[2][k] * dds2)
                                 + (xf->vs1[4][k] + xf->vs2[4][k] + xf->vsx[k]) * (xi_ule1 * dule1 + xi_ule2 * dule2);
            }

            if (b == xf->iblte[is] + 1) {
                /* redefine the coefficients for TTE, DTE, etc. */
                for (int k = 0; k < 3; k++) {
                    xf->vz[0][k] = xf->vs1[0][k] * cte_cte1;
                    xf->vz[1][k] = xf->vs1[0][k] * cte_tte1 + xf->vs1[1][k] * tte_tte1;
                    xf->vb[iv][0][k] = xf->vs1[0][k] * cte_cte2;
                    xf->vb[iv][1][k] = xf->vs1[0][k] * cte_tte2 + xf->vs1[1][k] * tte_tte2;
                }
            }

            /* turbulent intervals follow if currently at the transition interval */
            if (xf->tran) {
                xf->turb = 1;

                /* save the transition location */
                xf->itran[is] = b;
                xf->tforce[is] = xf->trforc;

                /* interpolate the airfoil geometry to find transition x/c (user output) */
                const double str = (is == 0) ? xf->sst - xf->xt : xf->sst + xf->xt;
                const double chx = xf->xte - xf->xle;
                const double chy = xf->yte - xf->yle;
                const double chsq = P2(chx) + P2(chy);
                const double xtr = seval(str, xf->x, xf->xp, xf->s, xf->n);
                const double ytr = seval(str, xf->y, xf->yp, xf->s, xf->n);
                xf->xoctr[is] = ((xtr - xf->xle) * chx + (ytr - xf->yle) * chy) / chsq;
            }

            xf->tran = 0;

            if (b == xf->iblte[is]) {
                /* "2" variables at the TE := wake correlations for the next station */
                xf->turb = 1;
                xf->wake = 1;
                blvar(xf, 3);
                blmid(xf, 3);
            }

            for (int js = 0; js < 2; js++) {
                for (int jbl = 1; jbl < xf->nbl[js]; jbl++) {
                    const int jv = xf->isys[js][jbl];
                    u1_m[jv] = u2_m[jv];
                    d1_m[jv] = d2_m[jv];
                }
            }

            u1_a = u2_a;
            d1_a = d2_a;

            due1 = due2;
            dds1 = dds2;

            /* BL variables for the next station */
            xf->com1 = xf->com2;
        }

        if (xfoil_print) {
            printf(" Side%2d%s transition at x/c = %7.4f%5d\n", is + 1, xf->tforce[is] ? " forced" : "  free ",
                   xf->xoctr[is], xf->itran[is] + 1);
        }
    }

    return result;
}

/* ------------------------------------------------------------------------- */
/* Part 5: update. */

/* a relaxation limit: accepted only if it reduces rlx to a positive value */
static double limit_rlx(double rlx, double candidate)
{
    return (candidate > 0.0 && candidate < rlx) ? candidate : rlx;
}

void update(xfoil_t *xf)
{
    /* new Ue (per BL station) and panel speeds, and their sensitivities wrt alpha or CL.
     * (Fortran places these scratch arrays on copies of VA/VB by EQUIVALENCE and copies
     * them back, overwriting parts of VA/VB that setbl rebuilds; clean-up C4d/C5d: work
     * arrays in the context, VA/VB left untouched.) */
    double (*unew)[IVX] = xf->scr.unew, (*u_ac)[IVX] = xf->scr.u_ac;
    double *qnew = xf->scr.qnew, *q_ac = xf->scr.q_ac;

    /* max allowable alpha and CL changes per iteration */
    const double dalmax = 0.5 * xf->dtor;
    const double dalmin = -(0.5 * xf->dtor);
    const double dclmax = 0.5;
    double dclmin = -0.5;
    if (xf->matyp != 1) dclmin = fmax_f(-0.5, -(0.9 * xf->cl));

    const double hstinv_u = xf->gamm1 * P2(xf->minf / xf->qinf) / (1.0 + 0.5 * xf->gamm1 * P2(xf->minf));

    /* new Ue distribution assuming no under-relaxation, and its sensitivity wrt alpha or Re */
    for (int is = 0; is < 2; is++) {
        for (int b = 1; b < xf->nbl[is]; b++) {
            const int i = xf->ipan[is][b];

            double dui = 0.0;
            double dui_ac = 0.0;
            for (int js = 0; js < 2; js++) {
                for (int jbl = 1; jbl < xf->nbl[js]; jbl++) {
                    const int j = xf->ipan[js][jbl];
                    const int jv = xf->isys[js][jbl];
                    const double ue_m = -(xf->vti[is][b] * xf->vti[js][jbl] * xf->dij[j][i]);
                    dui = dui + ue_m * (xf->mass[js][jbl] + xf->vdel[jv][0][2]);
                    dui_ac = dui_ac + ue_m * (-xf->vdel[jv][1][2]);
                }
            }

            /* UINV depends on "AC" only if "AC" is alpha */
            const double uinv_ac = xf->lalfa ? 0.0 : xf->uinv_a[is][b];

            unew[is][b] = xf->uinv[is][b] + dui;
            u_ac[is][b] = uinv_ac + dui_ac;
        }
    }

    /* new Qtan from the new Ue, with the appropriate sign change */
    for (int is = 0; is < 2; is++) {
        for (int b = 1; b <= xf->iblte[is]; b++) {
            const int i = xf->ipan[is][b];
            qnew[i] = xf->vti[is][b] * unew[is][b];
            q_ac[i] = xf->vti[is][b] * u_ac[is][b];
        }
    }

    /* new CL from this new Qtan */
    const double sa = sin(xf->alfa);
    const double ca = cos(xf->alfa);

    const double beta = sqrt(1.0 - P2(xf->minf));
    const double beta_msq = -(0.5 / beta);

    const double bfac = 0.5 * P2(xf->minf) / (1.0 + beta);
    const double bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    double clnew = 0.0, cl_a = 0.0, cl_ms = 0.0, cl_ac = 0.0;

    double cginc = 1.0 - P2(qnew[0] / xf->qinf);
    double cpg1 = cginc / (beta + bfac * cginc);
    double cpg1_ms = -(cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc));

    double cpi_q = -(2.0 * qnew[0] / P2(xf->qinf));
    double cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    double cpg1_ac = cpc_cpi * cpi_q * q_ac[0];

    for (int i = 0; i < xf->n; i++) {
        const int ip = (i == xf->n - 1) ? 0 : i + 1;

        cginc = 1.0 - P2(qnew[ip] / xf->qinf);
        const double cpg2 = cginc / (beta + bfac * cginc);
        const double cpg2_ms = -(cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc));

        cpi_q = -(2.0 * qnew[ip] / P2(xf->qinf));
        cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        const double cpg2_ac = cpc_cpi * cpi_q * q_ac[ip];

        const double dx = (xf->x[ip] - xf->x[i]) * ca + (xf->y[ip] - xf->y[i]) * sa;
        const double dx_a = -((xf->x[ip] - xf->x[i]) * sa) + (xf->y[ip] - xf->y[i]) * ca;

        const double ag = 0.5 * (cpg2 + cpg1);
        const double ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
        const double ag_ac = 0.5 * (cpg2_ac + cpg1_ac);

        clnew = clnew + dx * ag;
        cl_a = cl_a + dx_a * ag;
        cl_ms = cl_ms + dx * ag_ms;
        cl_ac = cl_ac + dx * ag_ac;

        cpg1 = cpg2;
        cpg1_ms = cpg2_ms;
        cpg1_ac = cpg2_ac;
    }

    /* under-relaxation factor; a limit may only reduce it to a positive value. (Fortran: with
     * Re(CL) (MATYP 2/3) and CL < 0, DCLMIN = max(-0.5, -0.9 CL) is positive and the DCLMIN limit
     * gave a negative or huge RLX; clean-up C6c. Normal cases are unchanged.) */
    xf->rlx = 1.0;

    double dac;
    if (xf->lalfa) {
        /* alpha prescribed, AC is CL: change in Re for the changing CL, since Re = Re(CL) */
        dac = (clnew - xf->cl) / (1.0 - cl_ac - cl_ms * 2.0 * xf->minf * xf->minf_cl);
        if (xf->rlx * dac > dclmax) xf->rlx = limit_rlx(xf->rlx, dclmax / dac);
        if (xf->rlx * dac < dclmin) xf->rlx = limit_rlx(xf->rlx, dclmin / dac);
    } else {
        /* CL prescribed, AC is alpha: change in alpha to drive CL to the prescribed value */
        dac = (clnew - xf->clspec) / (0.0 - cl_ac - cl_a);
        if (xf->rlx * dac > dalmax) xf->rlx = limit_rlx(xf->rlx, dalmax / dac);
        if (xf->rlx * dac < dalmin) xf->rlx = limit_rlx(xf->rlx, dalmin / dac);
    }

    xf->rmsbl = 0.0;
    xf->rmxbl = 0.0;

    const double dhi = 1.5;
    const double dlo = -0.5;

    /* changes in BL variables, and under-relaxation if needed */
    for (int is = 0; is < 2; is++) {
        for (int b = 1; b < xf->nbl[is]; b++) {
            const int iv = xf->isys[is][b];

            /* changes without under-relaxation */
            const double dctau = xf->vdel[iv][0][0] - dac * xf->vdel[iv][1][0];
            const double dthet = xf->vdel[iv][0][1] - dac * xf->vdel[iv][1][1];
            const double dmass = xf->vdel[iv][0][2] - dac * xf->vdel[iv][1][2];
            const double duedg = unew[is][b] + dac * u_ac[is][b] - xf->uedg[is][b];
            const double ddstr = (dmass - xf->dstr[is][b] * duedg) / xf->uedg[is][b];

            /* normalized changes */
            double dn1 = 0.0;
            if (b < xf->itran[is]) dn1 = dctau / 10.0;
            if (b >= xf->itran[is]) dn1 = dctau / xf->ctau[is][b];
            const double dn2 = dthet / xf->thet[is][b];
            const double dn3 = ddstr / xf->dstr[is][b];
            const double dn4 = fabs(duedg) / 0.25;

            /* accumulate for the rms change */
            xf->rmsbl = xf->rmsbl + P2(dn1) + P2(dn2) + P2(dn3) + P2(dn4);

            /* does Ctau need under-relaxation? */
            const double rdn1 = xf->rlx * dn1;
            if (fabs(dn1) > fabs(xf->rmxbl)) {
                xf->rmxbl = dn1;
                if (b < xf->itran[is]) xf->vmxbl[0] = 'n';
                if (b >= xf->itran[is]) xf->vmxbl[0] = 'C';
                xf->imxbl = b;
                xf->ismxbl = is;
            }
            if (rdn1 > dhi) xf->rlx = dhi / dn1;
            if (rdn1 < dlo) xf->rlx = dlo / dn1;

            /* Theta */
            const double rdn2 = xf->rlx * dn2;
            if (fabs(dn2) > fabs(xf->rmxbl)) {
                xf->rmxbl = dn2;
                xf->vmxbl[0] = 'T';
                xf->imxbl = b;
                xf->ismxbl = is;
            }
            if (rdn2 > dhi) xf->rlx = dhi / dn2;
            if (rdn2 < dlo) xf->rlx = dlo / dn2;

            /* Dstar */
            const double rdn3 = xf->rlx * dn3;
            if (fabs(dn3) > fabs(xf->rmxbl)) {
                xf->rmxbl = dn3;
                xf->vmxbl[0] = 'D';
                xf->imxbl = b;
                xf->ismxbl = is;
            }
            if (rdn3 > dhi) xf->rlx = dhi / dn3;
            if (rdn3 < dlo) xf->rlx = dlo / dn3;

            /* Ue */
            const double rdn4 = xf->rlx * dn4;
            if (fabs(dn4) > fabs(xf->rmxbl)) {
                xf->rmxbl = duedg;
                xf->vmxbl[0] = 'U';
                xf->imxbl = b;
                xf->ismxbl = is;
            }
            if (rdn4 > dhi) xf->rlx = dhi / dn4;
            if (rdn4 < dlo) xf->rlx = dlo / dn4;
        }
    }

    /* true rms change (NBL are counts) */
    xf->rmsbl = sqrt(xf->rmsbl / (4.0 * (double)(xf->nbl[0] + xf->nbl[1])));

    if (xf->lalfa) {
        /* under-relaxed change in Reynolds number from the change in lift */
        xf->cl = xf->cl + xf->rlx * dac;
    } else {
        /* under-relaxed change in alpha */
        xf->alfa = xf->alfa + xf->rlx * dac;
        xf->adeg = xf->alfa / xf->dtor;
    }

    /* update the BL variables with the under-relaxed changes */
    for (int is = 0; is < 2; is++) {
        for (int b = 1; b < xf->nbl[is]; b++) {
            const int iv = xf->isys[is][b];

            const double dctau = xf->vdel[iv][0][0] - dac * xf->vdel[iv][1][0];
            const double dthet = xf->vdel[iv][0][1] - dac * xf->vdel[iv][1][1];
            const double dmass = xf->vdel[iv][0][2] - dac * xf->vdel[iv][1][2];
            const double duedg = unew[is][b] + dac * u_ac[is][b] - xf->uedg[is][b];
            const double ddstr = (dmass - xf->dstr[is][b] * duedg) / xf->uedg[is][b];

            xf->ctau[is][b] = xf->ctau[is][b] + xf->rlx * dctau;
            xf->thet[is][b] = xf->thet[is][b] + xf->rlx * dthet;
            xf->dstr[is][b] = xf->dstr[is][b] + xf->rlx * ddstr;
            xf->uedg[is][b] = xf->uedg[is][b] + xf->rlx * duedg;

            double dswaki;
            if (b > xf->iblte[is]) {
                const int iw = b - xf->iblte[is];          /* Fortran IW = IBL - IBLTE (same value) */
                dswaki = xf->wgap[iw - 1];
            } else {
                dswaki = 0.0;
            }

            /* eliminate absurd transients */
            if (b >= xf->itran[is]) xf->ctau[is][b] = fmin_f(xf->ctau[is][b], 0.25);

            const double hklim = (b <= xf->iblte[is]) ? 1.02 : 1.00005;
            const double msq = P2(xf->uedg[is][b]) * hstinv_u / (xf->gamm1 * (1.0 - 0.5 * P2(xf->uedg[is][b]) * hstinv_u));
            double dsw = xf->dstr[is][b] - dswaki;
            dslim(&dsw, xf->thet[is][b], xf->uedg[is][b], msq, hklim);
            xf->dstr[is][b] = dsw + dswaki;

            /* new mass defect (nonlinear update) */
            xf->mass[is][b] = xf->dstr[is][b] * xf->uedg[is][b];
        }

        /* no "islands" of negative Ue (Fortran stations IBL = 3..IBLTE) */
        for (int b = 2; b <= xf->iblte[is]; b++) {
            if (xf->uedg[is][b - 1] > 0.0 && xf->uedg[is][b] <= 0.0) {
                xf->uedg[is][b] = xf->uedg[is][b - 1];
                xf->mass[is][b] = xf->dstr[is][b] * xf->uedg[is][b];
            }
        }
    }

    /* upper wake arrays := lower wake arrays (Fortran KBL = 1..NBL(2)-IBLTE(2)) */
    for (int kbl = 1; kbl <= xf->nbl[1] - (xf->iblte[1] + 1); kbl++) {
        xf->ctau[0][xf->iblte[0] + kbl] = xf->ctau[1][xf->iblte[1] + kbl];
        xf->thet[0][xf->iblte[0] + kbl] = xf->thet[1][xf->iblte[1] + kbl];
        xf->dstr[0][xf->iblte[0] + kbl] = xf->dstr[1][xf->iblte[1] + kbl];
        xf->uedg[0][xf->iblte[0] + kbl] = xf->uedg[1][xf->iblte[1] + kbl];
        xf->tau[0][xf->iblte[0] + kbl] = xf->tau[1][xf->iblte[1] + kbl];
        xf->dis[0][xf->iblte[0] + kbl] = xf->dis[1][xf->iblte[1] + kbl];
        xf->ctq[0][xf->iblte[0] + kbl] = xf->ctq[1][xf->iblte[1] + kbl];
        xf->delt[0][xf->iblte[0] + kbl] = xf->delt[1][xf->iblte[1] + kbl];
        xf->tstr[0][xf->iblte[0] + kbl] = xf->tstr[1][xf->iblte[1] + kbl];
    }

}
