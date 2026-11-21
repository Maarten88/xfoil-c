/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of s_xfoil.f90. max site in mrcl: standard order (objdump).
 */
#include "sxfoil.h"

#include <math.h>
#include <stdio.h>

#include "xstate.h"
#include "xfortran.h"

#define P2(v) __builtin_powi((v), 2)

void mrcl(xfoil_t *xf, double cls, double *m_cls, double *r_cls)
{
    const double cla = fmax_f(cls, 0.000001);

    if (xf->retyp < 1 || xf->retyp > 3) {
        if (xfoil_print) {
            printf(" MRCL:  Illegal Re(CL) dependence trigger.\n        Setting fixed Re.\n");
        }
        xf->retyp = 1;
    }
    if (xf->matyp < 1 || xf->matyp > 3) {
        if (xfoil_print) {
            printf(" MRCL:  Illegal Mach(CL) dependence trigger.\n        Setting fixed Mach.\n");
        }
        xf->matyp = 1;
    }

    if (xf->matyp == 1) {
        xf->minf = xf->minf1;
        *m_cls = 0.0;
    } else if (xf->matyp == 2) {
        xf->minf = xf->minf1 / sqrt(cla);
        *m_cls = -(0.5 * xf->minf / cla);
    } else if (xf->matyp == 3) {
        xf->minf = xf->minf1;
        *m_cls = 0.0;
    }

    if (xf->retyp == 1) {
        xf->reinf = xf->reinf1;
        *r_cls = 0.0;
    } else if (xf->retyp == 2) {
        xf->reinf = xf->reinf1 / sqrt(cla);
        *r_cls = -(0.5 * xf->reinf / cla);
    } else if (xf->retyp == 3) {
        xf->reinf = xf->reinf1 / cla;
        *r_cls = -(xf->reinf / cla);
    }

    if (xf->minf >= 0.99) {
        if (xfoil_print) {
            printf("\n MRCL: CL too low for chosen Mach(CL) dependence\n       Aritificially limiting Mach to  0.99\n");
        }
        xf->minf = 0.99;
        *m_cls = 0.0;
    }

    double rrat = 1.0;
    if (xf->reinf1 > 0.0) rrat = xf->reinf / xf->reinf1;

    if (rrat > 100.0) {
        if (xfoil_print) {
            printf("\n MRCL: CL too low for chosen Re(CL) dependence\n       Aritificially limiting Re to  %.17g\n", xf->reinf1 * 100.0);
        }
        xf->reinf = xf->reinf1 * 100.0;
        *r_cls = 0.0;
    }
}

void comset(xfoil_t *xf)
{
    /* Karman-Tsien parameter TKLAM */
    const double beta = sqrt(1.0 - P2(xf->minf));
    const double beta_msq = -(0.5 / beta);

    xf->tklam = P2(xf->minf) / P2(1.0 + beta);
    xf->tkl_msq = 1.0 / P2(1.0 + beta) - 2.0 * xf->tklam / (1.0 + beta) * beta_msq;

    /* sonic pressure coefficient and speed */
    if (xf->minf == 0.0) {
        xf->cpstar = -999.0;
        xf->qstar = 999.0;
    } else {
        xf->cpstar = 2.0 / (xf->gamma * P2(xf->minf))
                 * (pow((1.0 + 0.5 * xf->gamm1 * P2(xf->minf)) / (1.0 + 0.5 * xf->gamm1), xf->gamma / xf->gamm1) - 1.0);
        xf->qstar = xf->qinf / xf->minf * sqrt((1.0 + 0.5 * xf->gamm1 * P2(xf->minf)) / (1.0 + 0.5 * xf->gamm1));
    }
}

void cpcalc(int n, const double *q, double qinf, double minf, double *cp)
{
    const double beta = sqrt(1.0 - P2(minf));
    const double bfac = 0.5 * P2(minf) / (1.0 + beta);

    int denneg = 0;
    for (int i = 0; i < n; i++) {
        const double cpinc = 1.0 - P2(q[i] / qinf);
        const double den = beta + bfac * cpinc;
        cp[i] = cpinc / den;
        if (den <= 0.0) denneg = 1;
    }

    if (denneg && xfoil_print) {
        printf("\n CPCALC: Local speed too large. Compressibility corrections invalid.\n");
    }
}

void clcalc(int n, const double *x, const double *y, const double *gam, const double *gam_a,
            double alfa, double minf, double qinf, double xref, double yref,
            double *cl, double *cm, double *cdp, double *cl_alf, double *cl_msq)
{
    const double sa = sin(alfa);
    const double ca = cos(alfa);

    const double beta = sqrt(1.0 - P2(minf));
    const double beta_msq = -(0.5 / beta);

    const double bfac = 0.5 * P2(minf) / (1.0 + beta);
    const double bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    *cl = 0.0;
    *cm = 0.0;
    *cdp = 0.0;
    *cl_alf = 0.0;
    *cl_msq = 0.0;

    int i = 0;
    double cginc = 1.0 - P2(gam[i] / qinf);
    double cpg1 = cginc / (beta + bfac * cginc);
    double cpg1_msq = -(cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc));

    double cpi_gam = -(2.0 * gam[i] / P2(qinf));
    double cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    double cpg1_alf = cpc_cpi * cpi_gam * gam_a[i];

    for (i = 0; i < n; i++) {
        const int ip = (i == n - 1) ? 0 : i + 1;

        cginc = 1.0 - P2(gam[ip] / qinf);
        const double cpg2 = cginc / (beta + bfac * cginc);
        const double cpg2_msq = -(cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc));

        cpi_gam = -(2.0 * gam[ip] / P2(qinf));
        cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        const double cpg2_alf = cpc_cpi * cpi_gam * gam_a[ip];

        const double dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa;
        const double dy = (y[ip] - y[i]) * ca - (x[ip] - x[i]) * sa;
        const double dg = cpg2 - cpg1;

        const double ax = (0.5 * (x[ip] + x[i]) - xref) * ca + (0.5 * (y[ip] + y[i]) - yref) * sa;
        const double ay = (0.5 * (y[ip] + y[i]) - yref) * ca - (0.5 * (x[ip] + x[i]) - xref) * sa;
        const double ag = 0.5 * (cpg2 + cpg1);

        const double dx_alf = -((x[ip] - x[i]) * sa) + (y[ip] - y[i]) * ca;
        const double ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
        const double ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

        *cl = *cl + dx * ag;
        *cdp = *cdp - dy * ag;
        *cm = *cm - dx * (ag * ax + dg * dx / 12.0) - dy * (ag * ay + dg * dy / 12.0);

        *cl_alf = *cl_alf + dx * ag_alf + ag * dx_alf;
        *cl_msq = *cl_msq + dx * ag_msq;

        cpg1 = cpg2;
        cpg1_alf = cpg2_alf;
        cpg1_msq = cpg2_msq;
    }
}

void cdcalc(xfoil_t *xf)
{
    const double sa = sin(xf->alfa);
    const double ca = cos(xf->alfa);

    if (xf->lvisc && xf->lblini) {
        /* variables at the end of the wake (last station, 0-based nbl-1) */
        const int last = xf->nbl[1] - 1;
        const double thwake = xf->thet[1][last];
        const double urat = xf->uedg[1][last] / xf->qinf;
        const double uewake = xf->uedg[1][last] * (1.0 - xf->tklam) / (1.0 - xf->tklam * P2(urat));
        const double shwake = xf->dstr[1][last] / xf->thet[1][last];

        /* extrapolate the wake to downstream infinity (Squire-Young) */
        xf->cd = 2.0 * thwake * pow(uewake / xf->qinf, 0.5 * (5.0 + shwake));
    } else {
        xf->cd = 0.0;
    }

    /* friction drag coefficient (Fortran stations IBL = 3..IBLTE) */
    xf->cdf = 0.0;
    for (int is = 0; is < 2; is++) {
        for (int ibl = 2; ibl <= xf->iblte[is]; ibl++) {
            const int i = xf->ipan[is][ibl];
            const int im = xf->ipan[is][ibl - 1];
            const double dx = (xf->x[i] - xf->x[im]) * ca + (xf->y[i] - xf->y[im]) * sa;
            xf->cdf = xf->cdf + 0.5 * (xf->tau[is][ibl] + xf->tau[is][ibl - 1]) * dx * 2.0 / P2(xf->qinf);
        }
    }
}

void tecalc(xfoil_t *xf)
{
    /* TE base vector and TE bisector components */
    const double dxte = xf->x[0] - xf->x[xf->n - 1];
    const double dyte = xf->y[0] - xf->y[xf->n - 1];
    const double dxs = 0.5 * (-xf->xp[0] + xf->xp[xf->n - 1]);
    const double dys = 0.5 * (-xf->yp[0] + xf->yp[xf->n - 1]);

    /* normal and streamwise projected TE gap areas */
    xf->ante = dxs * dyte - dys * dxte;
    xf->aste = dxs * dxte + dys * dyte;

    /* total TE gap area */
    xf->dste = sqrt(P2(dxte) + P2(dyte));

    xf->sharp = xf->dste < 0.0001 * xf->chord;

    double scs, sds;
    if (xf->sharp) {
        scs = 1.0;
        sds = 0.0;
    } else {
        scs = xf->ante / xf->dste;
        sds = xf->aste / xf->dste;
    }

    /* TE panel source and vorticity strengths */
    xf->sigte = 0.5 * (xf->gam[0] - xf->gam[xf->n - 1]) * scs;
    xf->gamte = -(0.5 * (xf->gam[0] - xf->gam[xf->n - 1]) * sds);

}
