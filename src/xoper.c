/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of the API-reachable part of m_xoper.f90.
 */
#include "xoper.h"

#include <math.h>
#include <stdio.h>

#include "xstate.h"
#include "sxbl.h"
#include "sxfoil.h"
#include "sxoper.h"
#include "xbl.h"
#include "xfortran.h"
#include "xpanel.h"
#include "xsolve.h"

#define P2(v) __builtin_powi((v), 2)

void fcpmin(xfoil_t *xf)
{
    /* inviscid and viscous minimum Cp (Fortran globals CPMNI, CPMNV; only CPMN is used) */
    double cpmni = xf->cpi[0];
    double cpmnv = xf->cpv[0];

    for (int i = 1; i < xf->n + xf->nw; i++) {
        if (xf->cpi[i] < cpmni) {
            cpmni = xf->cpi[i];
        }
        if (xf->cpv[i] < cpmnv) {
            cpmnv = xf->cpv[i];
        }
    }

    if (xf->lvisc) {
        xf->cpmn = cpmnv;
    } else {
        xf->cpmn = cpmni;
    }
}

/* GAM, GAM_a, PSIO for the current alpha from the alpha = 0, 90 distributions. */
static void superimpose_gamu(xfoil_t *xf)
{
    xf->cosa = cos(xf->alfa);
    xf->sina = sin(xf->alfa);
    for (int i = 0; i < xf->n; i++) {
        xf->gam[i] = xf->cosa * xf->gamu[0][i] + xf->sina * xf->gamu[1][i];
        xf->gam_a[i] = -(xf->sina * xf->gamu[0][i]) + xf->cosa * xf->gamu[1][i];
    }
}

void specal(xfoil_t *xf)
{
    double minf_clm, reinf_clm;

    /* calculate surface vorticity distributions for alpha = 0, 90 degrees */
    if (!xf->lgamu || !xf->lqaij) ggcalc(xf);

    /* superimpose suitably weighted alpha = 0, 90 distributions */
    superimpose_gamu(xf);

    tecalc(xf);
    qiset(xf);

    /* set initial guess for the Newton variable CLM */
    double clm = 1.0;

    /* set corresponding M(CLM), Re(CLM) */
    mrcl(xf, clm, &minf_clm, &reinf_clm);
    comset(xf);

    /* set corresponding CL(M) */
    clcalc(xf->n, xf->x, xf->y, xf->gam, xf->gam_a, xf->alfa, xf->minf, xf->qinf, xf->xcmref, xf->ycmref, &xf->cl, &xf->cm, &xf->cdp, &xf->cl_alf, &xf->cl_msq);

    /* iterate on CLM */
    int converged = 0;
    for (int itcl = 1; itcl <= 20; itcl++) {
        const double msq_clm = 2.0 * xf->minf * minf_clm;
        const double dclm = (xf->cl - clm) / (1.0 - xf->cl_msq * msq_clm);

        const double clm1 = clm;
        xf->rlx = 1.0;

        /* under-relaxation loop to avoid driving M(CL) above 1 */
        for (int irlx = 1; irlx <= 12; irlx++) {
            clm = clm1 + xf->rlx * dclm;
            /* set new freestream Mach M(CLM) */
            mrcl(xf, clm, &minf_clm, &reinf_clm);
            /* if Mach is OK, go do next Newton iteration */
            if (xf->matyp == 1 || xf->minf == 0.0 || minf_clm != 0.0) break;
            xf->rlx = 0.5 * xf->rlx;
        }

        /* set new CL(M) */
        comset(xf);
        clcalc(xf->n, xf->x, xf->y, xf->gam, xf->gam_a, xf->alfa, xf->minf, xf->qinf, xf->xcmref, xf->ycmref, &xf->cl, &xf->cm, &xf->cdp, &xf->cl_alf, &xf->cl_msq);

        if (fabs(dclm) <= 1.0e-6) {
            converged = 1;
            break;
        }
    }
    if (!converged && xfoil_print) printf(" SPECAL:  Minf convergence failed\n");

    /* set final Mach, CL, Cp distributions, and hinge moment */
    mrcl(xf, xf->cl, &xf->minf_cl, &xf->reinf_cl);
    comset(xf);
    clcalc(xf->n, xf->x, xf->y, xf->gam, xf->gam_a, xf->alfa, xf->minf, xf->qinf, xf->xcmref, xf->ycmref, &xf->cl, &xf->cm, &xf->cdp, &xf->cl_alf, &xf->cl_msq);
    cpcalc(xf->n, xf->qinv, xf->qinf, xf->minf, xf->cpi);
    if (xf->lvisc) {
        cpcalc(xf->n + xf->nw, xf->qvis, xf->qinf, xf->minf, xf->cpv);
        cpcalc(xf->n + xf->nw, xf->qinv, xf->qinf, xf->minf, xf->cpi);
    } else {
        cpcalc(xf->n, xf->qinv, xf->qinf, xf->minf, xf->cpi);
    }
    if (xf->lflap) mhinge(xf);
}

void speccl(xfoil_t *xf)
{
    /* calculate surface vorticity distributions for alpha = 0, 90 degrees */
    if (!xf->lgamu || !xf->lqaij) ggcalc(xf);

    /* set freestream Mach from specified CL -- Mach will be held fixed */
    mrcl(xf, xf->clspec, &xf->minf_cl, &xf->reinf_cl);
    comset(xf);

    /* current alpha is the initial guess for Newton variable ALFA */
    superimpose_gamu(xf);

    /* get corresponding CL, CL_alpha, CL_Mach */
    clcalc(xf->n, xf->x, xf->y, xf->gam, xf->gam_a, xf->alfa, xf->minf, xf->qinf, xf->xcmref, xf->ycmref, &xf->cl, &xf->cm, &xf->cdp, &xf->cl_alf, &xf->cl_msq);

    /* Newton loop for alpha to get specified inviscid CL */
    int converged = 0;
    for (int ital = 1; ital <= 20; ital++) {
        const double dalfa = (xf->clspec - xf->cl) / xf->cl_alf;
        xf->rlx = 1.0;

        xf->alfa = xf->alfa + xf->rlx * dalfa;

        /* set new surface speed distribution */
        superimpose_gamu(xf);

        /* set new CL(alpha) */
        clcalc(xf->n, xf->x, xf->y, xf->gam, xf->gam_a, xf->alfa, xf->minf, xf->qinf, xf->xcmref, xf->ycmref, &xf->cl, &xf->cm, &xf->cdp, &xf->cl_alf, &xf->cl_msq);

        if (fabs(dalfa) <= 1.0e-6) {
            converged = 1;
            break;
        }
    }
    if (!converged && xfoil_print) printf(" SPECCL:  CL convergence failed\n");

    /* set final surface speed and Cp distributions */
    tecalc(xf);
    qiset(xf);
    if (xf->lvisc) {
        cpcalc(xf->n + xf->nw, xf->qvis, xf->qinf, xf->minf, xf->cpv);
        cpcalc(xf->n + xf->nw, xf->qinv, xf->qinf, xf->minf, xf->cpi);
    } else {
        cpcalc(xf->n, xf->qinv, xf->qinf, xf->minf, xf->cpi);
    }
    if (xf->lflap) mhinge(xf);
}

int viscal(xfoil_t *xf, int niter1)
{
    /* convergence tolerance */
    const double eps1 = 1.0e-4;
    int ok = 1;

    const int niter = niter1;   /* Fortran asks for the count if it is 0 (aski, not ported) */

    /* calculate wake trajectory from current inviscid solution if necessary */
    if (!xf->lwake) xywake(xf);

    /* set velocities on wake from airfoil vorticity for alpha=0, 90 */
    qwcalc(xf);

    /* set velocities on airfoil and wake for initial alpha */
    qiset(xf);

    if (!xf->lipan) {
        if (xf->lblini) gamqv(xf);
        /* locate stagnation point arc length position and panel index */
        stfind(xf);
        /* set BL position -> panel position pointers */
        iblpan(xf);
        /* calculate surface arc length array for current stagnation point location */
        xicalc(xf);
        /* set BL position -> system line pointers */
        iblsys(xf);
    }

    /* set inviscid BL edge velocity UINV from QINV */
    uicalc(xf);

    if (!xf->lblini) {
        /* set initial Ue from inviscid Ue */
        for (int b = 0; b < xf->nbl[0]; b++) xf->uedg[0][b] = xf->uinv[0][b];
        for (int b = 0; b < xf->nbl[1]; b++) xf->uedg[1][b] = xf->uinv[1][b];
    }

    if (xf->lvconv) {
        /* set correct CL if converged point exists */
        qvfue(xf);
        if (xf->lvisc) {
            cpcalc(xf->n + xf->nw, xf->qvis, xf->qinf, xf->minf, xf->cpv);
            cpcalc(xf->n + xf->nw, xf->qinv, xf->qinf, xf->minf, xf->cpi);
        } else {
            cpcalc(xf->n, xf->qinv, xf->qinf, xf->minf, xf->cpi);
        }
        gamqv(xf);
        clcalc(xf->n, xf->x, xf->y, xf->gam, xf->gam_a, xf->alfa, xf->minf, xf->qinf, xf->xcmref, xf->ycmref, &xf->cl, &xf->cm, &xf->cdp, &xf->cl_alf, &xf->cl_msq);
        cdcalc(xf);
    }

    /* set up source influence matrix if it doesn't exist */
    if (!xf->lwdij || !xf->ladij) qdcalc(xf);

    /* Newton iteration for entire BL solution */
    if (xfoil_print) {
        printf("\n");
        printf(" Solving BL system ...\n");
    }
    int converged = 0;
    for (int iter = 1; iter <= niter; iter++) {
        /* fill Newton system for BL variables */
        ok = setbl(xf);
        if (ABORT_ON_NAN && !ok) return ok;

        /* solve Newton system with custom solver */
        blsolv(xf);

        /* update BL variables */
        update(xf);

        if (xf->lalfa) {
            /* set new freestream Mach, Re from new CL */
            mrcl(xf, xf->cl, &xf->minf_cl, &xf->reinf_cl);
            comset(xf);
        } else {
            /* set new inviscid speeds QINV and UINV for new alpha */
            qiset(xf);
            uicalc(xf);
        }

        /* calculate edge velocities QVIS(.) from UEDG(..) */
        qvfue(xf);
        /* set GAM distribution from QVIS */
        gamqv(xf);
        /* relocate stagnation point */
        stmove(xf);
        /* set updated CL,CD */
        clcalc(xf->n, xf->x, xf->y, xf->gam, xf->gam_a, xf->alfa, xf->minf, xf->qinf, xf->xcmref, xf->ycmref, &xf->cl, &xf->cm, &xf->cdp, &xf->cl_alf, &xf->cl_msq);
        cdcalc(xf);

        /* display changes and test for convergence (plain translation of the messages) */
        if (xfoil_print) {
            printf("\n %3d   rms: %10.4e   max: %10.4e   %c at %4d%3d", iter, xf->rmsbl, xf->rmxbl, xf->vmxbl[0],
                   xf->imxbl + 1, xf->ismxbl + 1);
            if (xf->rlx < 1.0) printf("   RLX:%6.3f", xf->rlx);
            printf("\n");
            printf("       a =%7.3f      CL =%8.4f\n", xf->alfa / xf->dtor, xf->cl);
            printf("      Cm =%8.4f     CD =%9.5f   =>   CDf =%9.5f    CDp =%9.5f\n", xf->cm, xf->cd, xf->cdf, xf->cd - xf->cdf);
        }

        if (xf->rmsbl < eps1) {
            xf->lvconv = 1;
            xf->avisc = xf->alfa;
            xf->mvisc = xf->minf;
            converged = 1;
            break;
        }
    }
    if (!converged && xfoil_print) printf(" VISCAL:  Convergence failed\n");

    cpcalc(xf->n + xf->nw, xf->qinv, xf->qinf, xf->minf, xf->cpi);
    cpcalc(xf->n + xf->nw, xf->qvis, xf->qinf, xf->minf, xf->cpv);
    if (xf->lflap) mhinge(xf);

    if (xfoil_print) {
        /* separation/reattachment summary of side 1 (only printed in Fortran too) */
        const int is = 0;
        double hkmax = 0.0, hkm = 0.0, psep = 0.0, patt = 0.0;
        for (int b = 1; b <= xf->iblte[is]; b++) {
            const double hki = xf->dstr[is][b] / xf->thet[is][b];
            hkmax = fmax_f(hki, hkmax);
            if (hkm < 4.0 && hki >= 4.0) {
                const double hfrac = (4.0 - hkm) / (hki - hkm);
                const double pdefm = P2(xf->uedg[is][b - 1]) * xf->thet[is][b - 1];
                const double pdefi = P2(xf->uedg[is][b]) * xf->thet[is][b];
                psep = pdefm * (1.0 - hfrac) + pdefi * hfrac;
            }
            if (hkm > 4.0 && hki < 4.0) {
                const double hfrac = (4.0 - hkm) / (hki - hkm);
                const double pdefm = P2(xf->uedg[is][b - 1]) * xf->thet[is][b - 1];
                const double pdefi = P2(xf->uedg[is][b]) * xf->thet[is][b];
                patt = pdefm * (1.0 - hfrac) + pdefi * hfrac;
            }
            hkm = hki;
        }
        const double delp = patt - psep;
        printf(" %10.3f%10.4f%11.6f%11.6f%11.6f%11.6f%10.4f     #\n", xf->acrit[is], hkmax, xf->cd, 2.0 * psep,
               2.0 * patt, 2.0 * delp, xf->xoctr[is]);
    }
    return ok;
}
