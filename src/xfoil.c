/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of the API-reachable part of m_xfoil.f90.
 */
#include "xfoil.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "xstate.h"
#include "sxfoil.h"
#include "xbl.h"
#include "xerror.h"
#include "xfortran.h"
#include "xgeom.h"
#include "xnaca.h"
#include "xpanel.h"
#include "xspline.h"
#include "xuserio.h"

/* Fortran assignment of a literal to character(len): copy and blank-pad. */
#define FSET(dst, lit) fort_set_string((dst), sizeof(dst), (lit))
static void fort_set_string(char *dst, size_t len, const char *lit)
{
    size_t k = strlen(lit);
    if (k > len) k = len;
    memcpy(dst, lit, k);
    memset(dst + k, ' ', len - k);
}

static void init_impl(xfoil_t *xf)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */

    xf->pi = 4.0 * atan(1.0);
    xf->hopi = 0.50 / xf->pi;
    xf->qopi = 0.25 / xf->pi;
    xf->dtor = xf->pi / 180.0;

    /* default Cp/Cv (air) */
    xf->gamma = 1.4;
    xf->gamm1 = xf->gamma - 1.0;

    /* set unity freestream speed */
    xf->qinf = 1.0;

    /* initialize freestream Mach number to zero */
    xf->matyp = 1;
    xf->minf1 = 0.0;

    xf->alfa = 0.0;
    xf->cosa = 1.0;
    xf->sina = 0.0;

    for (int i = 0; i < IQX; i++) {
        xf->gamu[0][i] = 0.0;
        xf->gamu[1][i] = 0.0;
        xf->gam[i] = 0.0;
        xf->gam_a[i] = 0.0;
    }

    xf->cl = 0.0;
    xf->cm = 0.0;
    xf->cd = 0.0;

    xf->sigte = 0.0;
    xf->gamte = 0.0;

    for (int i = 0; i < IZX; i++) xf->sig[i] = 0.0;

    xf->awake = 0.0;
    xf->avisc = 0.0;

    xf->yimage = -10.0;
    xf->limage = 0;

    xf->lgamu = 0;
    xf->lvisc = 0;
    xf->lwake = 0;
    xf->lblini = 0;
    xf->lipan = 0;
    xf->lqaij = 0;
    xf->ladij = 0;
    xf->lwdij = 0;
    xf->lvconv = 0;
    xf->lbflap = 0;
    xf->lflap = 0;

    /* input airfoil will not be normalized */
    xf->lnorm = 0;

    /* airfoil will not be forced symmetric */
    xf->lgsym = 0;

    /* buffer and current airfoil flap hinge coordinates */
    xf->xbf = 0.0;
    xf->ybf = 0.0;
    xf->xof = 0.0;
    xf->yof = 0.0;

    /* default paneling parameters */
    xf->npan = 160;
    xf->cvpar = 1.0;
    xf->cterat = 0.15;
    xf->ctrrat = 0.2;

    /* default paneling refinement zone x/c endpoints */
    xf->xsref1 = 1.0;
    xf->xsref2 = 1.0;
    xf->xpref1 = 1.0;
    xf->xpref2 = 1.0;

    /* default Cm reference location */
    xf->xcmref = 0.25;
    xf->ycmref = 0.0;

    /* default viscous parameters */
    xf->retyp = 1;
    xf->reinf1 = 0.0;
    xf->acrit[0] = 9.0;
    xf->acrit[1] = 9.0;
    xf->xstrip[0] = 1.0;
    xf->xstrip[1] = 1.0;
    xf->xoctr[0] = 1.0;
    xf->xoctr[1] = 1.0;
    xf->waklen = 1.0;

    xf->idamp = 0;

    /* set BL calibration parameters */
    blpini(xf);

    /* Newton iteration limit */
    xf->itmax = 20;

    /* drop tolerance for BL system solver */
    xf->vaccel = 0.01;

    xf->nname = 32;
    FSET(xf->name, "                                ");


    /* set MINF, REINF, based on current CL-dependence */
    mrcl(xf, 1.0, &xf->minf_cl, &xf->reinf_cl);

    /* set various compressibility parameters from MINF */
    comset(xf);
}

/* ------------------------------------------------------------------------- */

void naca(xfoil_t *xf, int ides1)
{
    /* number of points per side */
    const int nside = IQX / 3;

    int ides;
    if (ides1 <= 0) {
        /* Fortran: call aski(...) (console prompt, not ported; user decision) */
        if (xfoil_print) printf(" This designation not implemented.\n");
        return;
    }
    ides = ides1;

    int itype = 0;
    if (ides <= 25099) itype = 5;
    if (ides <= 9999) itype = 4;

    if (itype == 0) {
        if (xfoil_print) printf(" This designation not implemented.\n");
        return;
    }

    if (itype == 4) naca4(ides, xf->w1, xf->w2, xf->w3, nside, xf->xb, xf->yb, &xf->nb, xf->name, (int)sizeof xf->name);
    if (itype == 5) naca5(&ides, xf->w1, xf->w2, xf->w3, nside, xf->xb, xf->yb, &xf->nb, xf->name, (int)sizeof xf->name);
    strip(xf->name, (int)sizeof xf->name, &xf->nname);

    /* see if routines didn't recognize designator */
    if (ides == 0) return;


    xf->xbf = 0.0;
    xf->ybf = 0.0;
    xf->lbflap = 0;

    scalc(xf->xb, xf->yb, xf->sb, xf->nb);
    segspl(xf->xb, xf->xbp, xf->sb, xf->nb);
    segspl(xf->yb, xf->ybp, xf->sb, xf->nb);

    geopar(xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb, xf->w1, &xf->sble, &xf->chordb, &xf->areab, &xf->radble, &xf->angbte, &xf->ei11ba, &xf->ei22ba,
           &xf->apx1ba, &xf->apx2ba, &xf->ei11bt, &xf->ei22bt, &xf->apx1bt, &xf->apx2bt, &xf->thickb, &xf->cambrb);

    if (xfoil_print) printf("\n Buffer airfoil set using%4d points\n", xf->nb);

    /* set paneling */
    pangen(xf, 1);
}

#define P2(v) __builtin_powi((v), 2)

void pangen(xfoil_t *xf, int shopar)
{
    if (xf->nb < 2) {
        if (xfoil_print) printf(" PANGEN: Buffer airfoil not available.\n");
        xf->n = 0;
        return;
    }

    /* number of temporary nodes for the panel distribution calculation exceeds the
     * specified panel number by a factor of IPFAC (Fortran sets 3, then 5) */
    const int ipfac = 5;

    /* number of airfoil panel points */
    xf->n = xf->npan;

    /* set arc length spline parameter */
    scalc(xf->xb, xf->yb, xf->sb, xf->nb);

    /* spline raw airfoil coordinates */
    segspl(xf->xb, xf->xbp, xf->sb, xf->nb);
    segspl(xf->yb, xf->ybp, xf->sb, xf->nb);

    /* normalizing length (~ chord) */
    const double sbref = 0.5 * (xf->sb[xf->nb - 1] - xf->sb[0]);

    /* set up curvature array */
    for (int i = 0; i < xf->nb; i++) xf->w5[i] = fabs(curv(xf->sb[i], xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb)) * sbref;

    /* locate LE point arc length value and the normalized curvature there */
    lefind(&xf->sble, xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb);
    const double cvle = fabs(curv(xf->sble, xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb)) * sbref;

    /* check for doubled point (sharp corner) at LE: ible = first of the two nodes, -1 = none */
    int ible = -1;
    for (int i = 0; i < xf->nb - 1; i++) {
        if (xf->sble == xf->sb[i] && xf->sble == xf->sb[i + 1]) {
            ible = i;
            if (xfoil_print) {
                printf("\n");
                printf(" Sharp leading edge\n");
            }
            break;
        }
    }

    /* set LE, TE points */
    const double xble = seval(xf->sble, xf->xb, xf->xbp, xf->sb, xf->nb);
    const double yble = seval(xf->sble, xf->yb, xf->ybp, xf->sb, xf->nb);
    const double xbte = 0.5 * (xf->xb[0] + xf->xb[xf->nb - 1]);
    const double ybte = 0.5 * (xf->yb[0] + xf->yb[xf->nb - 1]);
    const double chbsq = P2(xbte - xble) + P2(ybte - yble);

    /* set average curvature over 2*NK+1 points within Rcurv of LE point */
    const int nk = 3;
    double cvsum = 0.0;
    for (int k = -nk; k <= nk; k++) {
        const double frac = (double)k / (double)nk;
        const double sbk = xf->sble + frac * sbref / fmax_f(cvle, 20.0);
        const double cvk = fabs(curv(sbk, xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb)) * sbref;
        cvsum = cvsum + cvk;
    }
    double cvavg = cvsum / (double)(2 * nk + 1);

    /* dummy curvature for sharp LE */
    if (ible >= 0) cvavg = 10.0;

    /* set curvature attraction coefficient actually used */
    const double cc = 6.0 * xf->cvpar;

    /* set artificial curvature at TE to bunch panels there */
    const double cvte = cvavg * xf->cterat;
    xf->w5[0] = cvte;
    xf->w5[xf->nb - 1] = cvte;

    /* smooth curvature array for smoother panel size distribution: smoothing length =
     * 1 / averaged LE curvature, but no more than 5% of chord and no less than 1/4
     * average panel spacing */
    const double smool = fmax_f(1.0 / fmax_f(cvavg, 20.0), 0.25 / (double)(xf->npan / 2));

    const double smoosq = P2(smool * sbref);

    /* set up tri-diagonal system for smoothed curvatures */
    xf->w2[0] = 1.0;
    xf->w3[0] = 0.0;
    for (int i = 1; i < xf->nb - 1; i++) {
        const double dsm = xf->sb[i] - xf->sb[i - 1];
        const double dsp = xf->sb[i + 1] - xf->sb[i];
        const double dso = 0.5 * (xf->sb[i + 1] - xf->sb[i - 1]);

        if (dsm == 0.0 || dsp == 0.0) {
            /* leave curvature at corner point unchanged */
            xf->w1[i] = 0.0;
            xf->w2[i] = 1.0;
            xf->w3[i] = 0.0;
        } else {
            xf->w1[i] = smoosq * (-1.0 / dsm) / dso;
            xf->w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            xf->w3[i] = smoosq * (-1.0 / dsp) / dso;
        }
    }

    xf->w1[xf->nb - 1] = 0.0;
    xf->w2[xf->nb - 1] = 1.0;

    /* fix curvature at LE point by modifying equations adjacent to LE */
    for (int i = 1; i < xf->nb - 1; i++) {
        if (xf->sb[i] == xf->sble || i == ible || i == ible + 1) {
            /* if node falls right on LE point, fix curvature there */
            xf->w1[i] = 0.0;
            xf->w2[i] = 1.0;
            xf->w3[i] = 0.0;
            xf->w5[i] = cvle;
        } else if (xf->sb[i - 1] < xf->sble && xf->sb[i] > xf->sble) {
            /* modify equation at node just before LE point; not at i = 1, where that node is
             * the TE end with its fixed equation (Fortran read SB(0) there; clean-up C6a) */
            double dsm, dsp, dso;
            if (i >= 2) {
                dsm = xf->sb[i - 1] - xf->sb[i - 2];
                dsp = xf->sble - xf->sb[i - 1];
                dso = 0.5 * (xf->sble - xf->sb[i - 2]);

                xf->w1[i - 1] = smoosq * (-1.0 / dsm) / dso;
                xf->w2[i - 1] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
                xf->w3[i - 1] = 0.0;
                xf->w5[i - 1] = xf->w5[i - 1] + smoosq * cvle / (dsp * dso);
            }

            /* modify equation at node just after LE point */
            dsm = xf->sb[i] - xf->sble;
            dsp = xf->sb[i + 1] - xf->sb[i];
            dso = 0.5 * (xf->sb[i + 1] - xf->sble);
            xf->w1[i] = 0.0;
            xf->w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            xf->w3[i] = smoosq * (-1.0 / dsp) / dso;
            xf->w5[i] = xf->w5[i] + smoosq * cvle / (dsm * dso);

            break;
        }
    }

    /* set artificial curvature at bunching points and fix it there */
    for (int i = 1; i < xf->nb - 1; i++) {
        /* chord-based x/c coordinate */
        const double xoc = ((xf->xb[i] - xble) * (xbte - xble) + (xf->yb[i] - yble) * (ybte - yble)) / chbsq;

        if (xf->sb[i] < xf->sble) {
            /* check if top side point is in refinement area */
            if (xoc > xf->xsref1 && xoc < xf->xsref2) {
                xf->w1[i] = 0.0;
                xf->w2[i] = 1.0;
                xf->w3[i] = 0.0;
                xf->w5[i] = cvle * xf->ctrrat;
            }
        } else if (xoc > xf->xpref1 && xoc < xf->xpref2) {
            /* bottom side point is in refinement area */
            xf->w1[i] = 0.0;
            xf->w2[i] = 1.0;
            xf->w3[i] = 0.0;
            xf->w5[i] = cvle * xf->ctrrat;
        }
    }

    /* solve for smoothed curvature array W5 */
    if (ible < 0) {
        trisol(xf->w2, xf->w1, xf->w3, xf->w5, xf->nb);
    } else {
        /* separately on both sides of the LE corner: nodes 0..ible and ible+1..nb-1 */
        trisol(xf->w2, xf->w1, xf->w3, xf->w5, ible + 1);
        const int i = ible + 1;
        trisol(&xf->w2[i], &xf->w1[i], &xf->w3[i], &xf->w5[i], xf->nb - i);
    }

    /* find max curvature */
    double cvmax = 0.0;
    for (int i = 0; i < xf->nb; i++) cvmax = fmax_f(cvmax, fabs(xf->w5[i]));

    /* normalize curvature array */
    for (int i = 0; i < xf->nb; i++) xf->w5[i] = xf->w5[i] / cvmax;

    /* spline curvature array */
    segspl(xf->w5, xf->w6, xf->sb, xf->nb);

    /* initial guess for node positions uniform in s; more nodes than specified
     * (by factor of IPFAC) are temporarily used for more reliable convergence */
    const int nn = ipfac * (xf->n - 1) + 1;

    /* ratio of lengths of panel at TE to one away from the TE */
    const double rdste = 0.667;
    const double rtf = (rdste - 1.0) * 2.0 + 1.0;

    int nn1 = 0;     /* number of nodes up to and including a sharp LE */
    if (ible < 0) {
        const double dsavg = (xf->sb[xf->nb - 1] - xf->sb[0]) / ((double)(nn - 3) + 2.0 * rtf);
        xf->snew[0] = xf->sb[0];
        for (int i = 1; i < nn - 1; i++) xf->snew[i] = xf->sb[0] + dsavg * ((double)(i - 1) + rtf);
        xf->snew[nn - 1] = xf->sb[xf->nb - 1];
    } else {
        const int nfrac1 = (xf->n * (ible + 1)) / xf->nb;

        nn1 = ipfac * (nfrac1 - 1) + 1;
        const double dsavg1 = (xf->sble - xf->sb[0]) / ((double)(nn1 - 2) + rtf);
        xf->snew[0] = xf->sb[0];
        for (int i = 1; i < nn1; i++) xf->snew[i] = xf->sb[0] + dsavg1 * ((double)(i - 1) + rtf);

        const int nn2 = nn - nn1 + 1;
        const double dsavg2 = (xf->sb[xf->nb - 1] - xf->sble) / ((double)(nn2 - 2) + rtf);
        for (int i = 1; i < nn2 - 1; i++) xf->snew[i + nn1 - 1] = xf->sble + dsavg2 * ((double)(i - 1) + rtf);
        xf->snew[nn - 1] = xf->sb[xf->nb - 1];
    }

    /* Newton iteration loop for new node positions */
    int converged = 0;
    for (int iter = 0; iter < 20; iter++) {
        /* set up tri-diagonal system for node position deltas */
        double cv1 = seval(xf->snew[0], xf->w5, xf->w6, xf->sb, xf->nb);
        double cv2 = seval(xf->snew[1], xf->w5, xf->w6, xf->sb, xf->nb);
        double cvs1 = deval(xf->snew[0], xf->w5, xf->w6, xf->sb, xf->nb);
        double cvs2 = deval(xf->snew[1], xf->w5, xf->w6, xf->sb, xf->nb);

        double cavm = sqrt(P2(cv1) + P2(cv2));
        double cavm_s1, cavm_s2;
        if (cavm == 0.0) {
            cavm_s1 = 0.0;
            cavm_s2 = 0.0;
        } else {
            cavm_s1 = cvs1 * cv1 / cavm;
            cavm_s2 = cvs2 * cv2 / cavm;
        }

        for (int i = 1; i < nn - 1; i++) {
            const double dsm = xf->snew[i] - xf->snew[i - 1];
            const double dsp = xf->snew[i] - xf->snew[i + 1];
            const double cv3 = seval(xf->snew[i + 1], xf->w5, xf->w6, xf->sb, xf->nb);
            const double cvs3 = deval(xf->snew[i + 1], xf->w5, xf->w6, xf->sb, xf->nb);

            const double cavp = sqrt(P2(cv3) + P2(cv2));
            double cavp_s2, cavp_s3;
            if (cavp == 0.0) {
                cavp_s2 = 0.0;
                cavp_s3 = 0.0;
            } else {
                cavp_s2 = cvs2 * cv2 / cavp;
                cavp_s3 = cvs3 * cv3 / cavp;
            }

            const double fm = cc * cavm + 1.0;
            const double fp = cc * cavp + 1.0;

            const double rez = dsp * fp + dsm * fm;

            /* lower, main, and upper diagonals */
            xf->w1[i] = -fm + cc * dsm * cavm_s1;
            xf->w2[i] = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2);
            xf->w3[i] = -fp + cc * dsp * cavp_s3;

            /* residual, requiring that (1 + C*curv)*deltaS is equal on both sides of node i */
            xf->w4[i] = -rez;

            cv2 = cv3;
            cvs2 = cvs3;
            cavm = cavp;
            cavm_s1 = cavp_s2;
            cavm_s2 = cavp_s3;
        }

        /* fix endpoints (at TE) */
        xf->w2[0] = 1.0;
        xf->w3[0] = 0.0;
        xf->w4[0] = 0.0;
        xf->w1[nn - 1] = 0.0;
        xf->w2[nn - 1] = 1.0;
        xf->w4[nn - 1] = 0.0;

        if (rtf != 1.0) {
            /* fudge equations adjacent to TE to get TE panel length ratio RTF */
            int i = 1;
            xf->w4[i] = -((xf->snew[i] - xf->snew[i - 1]) + rtf * (xf->snew[i] - xf->snew[i + 1]));
            xf->w1[i] = -1.0;
            xf->w2[i] = 1.0 + rtf;
            xf->w3[i] = -rtf;

            i = nn - 2;
            xf->w4[i] = -((xf->snew[i] - xf->snew[i + 1]) + rtf * (xf->snew[i] - xf->snew[i - 1]));
            xf->w3[i] = -1.0;
            xf->w2[i] = 1.0 + rtf;
            xf->w1[i] = -rtf;
        }

        /* fix sharp LE point */
        if (ible >= 0) {
            const int i = nn1 - 1;
            xf->w1[i] = 0.0;
            xf->w2[i] = 1.0;
            xf->w3[i] = 0.0;
            xf->w4[i] = xf->sble - xf->snew[i];
        }

        /* solve for changes W4 in node position arc length values */
        trisol(xf->w2, xf->w1, xf->w3, xf->w4, nn);

        /* find under-relaxation factor to keep nodes from changing order */
        xf->rlx = 1.0;
        double dmax = 0.0;
        for (int i = 0; i < nn - 1; i++) {
            const double ds = xf->snew[i + 1] - xf->snew[i];
            const double dds = xf->w4[i + 1] - xf->w4[i];
            const double dsrat = 1.0 + xf->rlx * dds / ds;
            if (dsrat > 4.0) xf->rlx = (4.0 - 1.0) * ds / dds;
            if (dsrat < 0.2) xf->rlx = (0.2 - 1.0) * ds / dds;
            dmax = fmax_f(fabs(xf->w4[i]), dmax);
        }

        /* update node position */
        for (int i = 1; i < nn - 1; i++) xf->snew[i] = xf->snew[i] + xf->rlx * xf->w4[i];

        if (fabs(dmax) < 1.e-3) {
            converged = 1;
            break;
        }
    }
    if (!converged && xfoil_print) printf(" Paneling convergence failed.  Continuing anyway...\n");

    /* set new panel node coordinates */
    for (int i = 0; i < xf->n; i++) {
        const int ind = ipfac * i;
        xf->s[i] = xf->snew[ind];
        xf->x[i] = seval(xf->snew[ind], xf->xb, xf->xbp, xf->sb, xf->nb);
        xf->y[i] = seval(xf->snew[ind], xf->yb, xf->ybp, xf->sb, xf->nb);
    }

    /* go over buffer airfoil again, checking for corners (double points) */
    for (int ib = 0; ib < xf->nb - 1; ib++) {
        if (xf->sb[ib] == xf->sb[ib + 1]) {
            /* found one! */
            const double xbcorn = xf->xb[ib];
            const double ybcorn = xf->yb[ib];
            const double sbcorn = xf->sb[ib];

            /* find current-airfoil panel which contains corner (Fortran trip count fixed at entry) */
            const int nloop = xf->n;
            for (int i = 0; i < nloop; i++) {
                /* keep stepping until first node past corner */
                if (xf->s[i] > sbcorn) {
                    /* move remainder of panel nodes to make room for additional node */
                    for (int j = xf->n; j > i; j--) {
                        xf->x[j] = xf->x[j - 1];
                        xf->y[j] = xf->y[j - 1];
                        xf->s[j] = xf->s[j - 1];
                    }
                    xf->n = xf->n + 1;

                    if (xf->n > IQX - 1) xfoil_stop("PANEL: Too many panels. Increase IQX in XFOIL.INC");

                    xf->x[i] = xbcorn;
                    xf->y[i] = ybcorn;
                    xf->s[i] = sbcorn;

                    /* shift nodes adjacent to corner to keep panel sizes comparable */
                    if (i >= 2) {
                        xf->s[i - 1] = 0.5 * (xf->s[i] + xf->s[i - 2]);
                        xf->x[i - 1] = seval(xf->s[i - 1], xf->xb, xf->xbp, xf->sb, xf->nb);
                        xf->y[i - 1] = seval(xf->s[i - 1], xf->yb, xf->ybp, xf->sb, xf->nb);
                    }

                    if (i + 2 < xf->n) {
                        xf->s[i + 1] = 0.5 * (xf->s[i] + xf->s[i + 2]);
                        xf->x[i + 1] = seval(xf->s[i + 1], xf->xb, xf->xbp, xf->sb, xf->nb);
                        xf->y[i + 1] = seval(xf->s[i + 1], xf->yb, xf->ybp, xf->sb, xf->nb);
                    }

                    /* go on to next input geometry point to check for corner */
                    break;
                }
            }
        }
    }

    scalc(xf->x, xf->y, xf->s, xf->n);
    segspl(xf->x, xf->xp, xf->s, xf->n);
    segspl(xf->y, xf->yp, xf->s, xf->n);
    lefind(&xf->sle, xf->x, xf->xp, xf->y, xf->yp, xf->s, xf->n);

    xf->xle = seval(xf->sle, xf->x, xf->xp, xf->s, xf->n);
    xf->yle = seval(xf->sle, xf->y, xf->yp, xf->s, xf->n);
    xf->xte = 0.5 * (xf->x[0] + xf->x[xf->n - 1]);
    xf->yte = 0.5 * (xf->y[0] + xf->y[xf->n - 1]);
    xf->chord = sqrt(P2(xf->xte - xf->xle) + P2(xf->yte - xf->yle));

    /* Fortran computes panel size ratios dsmin/dsmax here only for a commented-out
     * print; dropped (no effect on results). */

    /* set various flags for new airfoil */
    xf->lgamu = 0;
    xf->lwake = 0;
    xf->lqaij = 0;
    xf->ladij = 0;
    xf->lwdij = 0;
    xf->lipan = 0;
    xf->lblini = 0;
    xf->lvconv = 0;

    if (xf->lbflap) {
        xf->xof = xf->xbf;
        xf->yof = xf->ybf;
        xf->lflap = 1;
    }

    /* determine if TE is blunt or sharp, calculate TE geometry parameters */
    tecalc(xf);

    /* calculate normal vectors */
    ncalc(xf->x, xf->y, xf->s, xf->n, xf->nx, xf->ny);

    /* calculate panel angles for panel routines */
    apcalc(xf);

    if (xfoil_print) {
        if (xf->sharp) {
            printf("\n Sharp trailing edge\n");
        } else {
            const double gap = sqrt(P2(xf->x[0] - xf->x[xf->n - 1]) + P2(xf->y[0] - xf->y[xf->n - 1]));
            printf("\n Blunt trailing edge.  Gap =%9.5f\n", gap);
        }
        if (shopar) {
            printf("\n Paneling parameters used...\n");
            printf("   Number of panel nodes      %4d\n", xf->npan);
            printf("   Panel bunching parameter   %6.3f\n", xf->cvpar);
            printf("   TE/LE panel density ratio  %6.3f\n", xf->cterat);
            printf("   Refined-area/LE panel density ratio   %6.3f\n", xf->ctrrat);
            printf("   Top    side refined area x/c limits %6.3f%6.3f\n", xf->xsref1, xf->xsref2);
            printf("   Bottom side refined area x/c limits %6.3f%6.3f\n", xf->xpref1, xf->xpref2);
        }
    }
}

/* ---- public entry points: run inside the error guard (xerror.h) ---- */

void init(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(init_impl(xf));
}
