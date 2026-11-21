/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of the API-reachable routines of m_xpanel.f90 (part 1: apcalc, ncalc,
 *  psilin, pswlin). Expression order follows the Fortran. Fortran unary minus
 *  binds weaker than * and /: `-a*b` is -(a*b), written so in C (it differs
 *  from (-a)*b for NaN, whose sign bit only negation flips).
 */
#include "xpanel.h"

#include <math.h>

#include <stdio.h>

#include "xstate.h"
#include "sxbl.h"
#include "xerror.h"
#include "xfortran.h"
#include "xsolve.h"
#include "xspline.h"
#include "xutils.h"

#define SQ(v) __builtin_powi((v), 2)    /* Fortran (v)**2 */

void apcalc(xfoil_t *xf)
{
    /* angles of airfoil panels */
    for (int i = 0; i < xf->n - 1; i++) {
        const double sx = xf->x[i + 1] - xf->x[i];
        const double sy = xf->y[i + 1] - xf->y[i];
        if (sx == 0.0 && sy == 0.0)
            xf->apanel[i] = atan2(-xf->ny[i], -xf->nx[i]);
        else
            xf->apanel[i] = atan2(sx, -sy);
    }

    /* TE panel */
    const int i = xf->n - 1;
    const int ip = 0;
    if (xf->sharp) {
        xf->apanel[i] = xf->pi;
    } else {
        const double sx = xf->x[ip] - xf->x[i];
        const double sy = xf->y[ip] - xf->y[i];
        xf->apanel[i] = atan2(-sx, sy) + xf->pi;
    }
}

void ncalc(const double *x, const double *y, const double *s, int n, double *xn, double *yn)
{
    if (n <= 1) return;

    segspl(x, xn, s, n);
    segspl(y, yn, s, n);
    for (int i = 0; i < n; i++) {
        const double sx = yn[i];
        const double sy = -xn[i];
        const double smod = sqrt(sx * sx + sy * sy);
        if (smod == 0.0) {
            xn[i] = -1.0;
            yn[i] = 0.0;
        } else {
            xn[i] = sx / smod;
            yn[i] = sy / smod;
        }
    }

    /* average normal vectors at corner points */
    for (int i = 0; i < n - 1; i++) {
        if (s[i] == s[i + 1]) {
            const double sx = 0.5 * (xn[i] + xn[i + 1]);
            const double sy = 0.5 * (yn[i] + yn[i + 1]);
            const double smod = sqrt(sx * sx + sy * sy);
            if (smod == 0.0) {
                xn[i] = -1.0;
                yn[i] = 0.0;
                xn[i + 1] = -1.0;
                yn[i + 1] = 0.0;
            } else {
                xn[i] = sx / smod;
                yn[i] = sy / smod;
                xn[i + 1] = sx / smod;
                yn[i + 1] = sy / smod;
            }
        }
    }
}

/* Linear source contributions of the 1-0 and 0-2 half-panels of panel jo-jp
 * (shared by psilin's real and image loops and by pswlin). tsign is the sign
 * of the (0.5 - 0.5*sgn)*pi branch shift (+1 in psilin, -1 in pswlin). */
static void source_half_panels(xfoil_t *xf, int jm, int jo, int jp, int jq, double x1, double x2, double yy,
                               double rs1, double rs2, double g1, double g2, double t1, double t2,
                               double apan, double sgn, int tsign, double x1i, double x2i, double yyi,
                               double dsio, double *psi, double *psi_ni)
{
    /* midpoint quantities */
    const double x0 = 0.5 * (x1 + x2);
    const double rs0 = x0 * x0 + yy * yy;
    const double g0 = log(rs0);
    const double shift = (0.5 - 0.5 * sgn) * xf->pi;
    const double t0 = (tsign > 0) ? atan2(sgn * x0, sgn * yy) + shift : atan2(sgn * x0, sgn * yy) - shift;

    /* source contribution to Psi for the 1-0 half-panel */
    double dxinv = 1.0 / (x1 - x0);
    double psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
    double pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;

    const double psx1 = -(t1 - apan);
    double psx0 = t0 - apan;
    double psyy = 0.5 * (g1 - g0);

    const double pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
    double pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
    double pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

    const double dsm = sqrt(SQ(xf->x[jp] - xf->x[jm]) + SQ(xf->y[jp] - xf->y[jm]));
    const double dsim = 1.0 / dsm;

    double ssum = (xf->sig[jp] - xf->sig[jo]) * dsio + (xf->sig[jp] - xf->sig[jm]) * dsim;
    double sdif = (xf->sig[jp] - xf->sig[jo]) * dsio - (xf->sig[jp] - xf->sig[jm]) * dsim;

    *psi = *psi + xf->qopi * (psum * ssum + pdif * sdif);

    /* dPsi/dm */
    xf->dzdm[jm] = xf->dzdm[jm] + xf->qopi * (-(psum * dsim) + pdif * dsim);
    xf->dzdm[jo] = xf->dzdm[jo] + xf->qopi * (-(psum * dsio) - pdif * dsio);
    xf->dzdm[jp] = xf->dzdm[jp] + xf->qopi * (psum * (dsio + dsim) + pdif * (dsio - dsim));

    /* dPsi/dni */
    double psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
    double pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
    *psi_ni = *psi_ni + xf->qopi * (psni * ssum + pdni * sdif);

    xf->dqdm[jm] = xf->dqdm[jm] + xf->qopi * (-(psni * dsim) + pdni * dsim);
    xf->dqdm[jo] = xf->dqdm[jo] + xf->qopi * (-(psni * dsio) - pdni * dsio);
    xf->dqdm[jp] = xf->dqdm[jp] + xf->qopi * (psni * (dsio + dsim) + pdni * (dsio - dsim));

    /* source contribution to Psi for the 0-2 half-panel */
    dxinv = 1.0 / (x0 - x2);
    psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
    pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;

    psx0 = -(t0 - apan);
    const double psx2 = t2 - apan;
    psyy = 0.5 * (g0 - g2);

    pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
    const double pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
    pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

    const double dsp = sqrt(SQ(xf->x[jq] - xf->x[jo]) + SQ(xf->y[jq] - xf->y[jo]));
    const double dsip = 1.0 / dsp;

    ssum = (xf->sig[jq] - xf->sig[jo]) * dsip + (xf->sig[jp] - xf->sig[jo]) * dsio;
    sdif = (xf->sig[jq] - xf->sig[jo]) * dsip - (xf->sig[jp] - xf->sig[jo]) * dsio;

    *psi = *psi + xf->qopi * (psum * ssum + pdif * sdif);

    /* dPsi/dm */
    xf->dzdm[jo] = xf->dzdm[jo] + xf->qopi * (-(psum * (dsip + dsio)) - pdif * (dsip - dsio));
    xf->dzdm[jp] = xf->dzdm[jp] + xf->qopi * (psum * dsio - pdif * dsio);
    xf->dzdm[jq] = xf->dzdm[jq] + xf->qopi * (psum * dsip + pdif * dsip);

    /* dPsi/dni */
    psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
    pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
    *psi_ni = *psi_ni + xf->qopi * (psni * ssum + pdni * sdif);

    xf->dqdm[jo] = xf->dqdm[jo] + xf->qopi * (-(psni * (dsip + dsio)) - pdni * (dsip - dsio));
    xf->dqdm[jp] = xf->dqdm[jp] + xf->qopi * (psni * dsio - pdni * dsio);
    xf->dqdm[jq] = xf->dqdm[jq] + xf->qopi * (psni * dsip + pdni * dsip);
}

/* Geometric sensitivities of x1, x2, yy to the normal motion of nodes jo, jp */
struct geo { double x1o, x1p, x2o, x2p, yyo, yyp; };

static struct geo geo_derivs(xfoil_t *xf, int jo, int jp, double rx1, double ry1, double rx2, double ry2,
                             double x1, double x2, double sx, double sy, double dsio)
{
    const double nxo = xf->nx[jo], nyo = xf->ny[jo], nxp = xf->nx[jp], nyp = xf->ny[jp];
    struct geo g;
    g.x1o = -(((rx1 - x1 * sx) * nxo + (ry1 - x1 * sy) * nyo) * dsio) - (sx * nxo + sy * nyo);
    g.x1p = ((rx1 - x1 * sx) * nxp + (ry1 - x1 * sy) * nyp) * dsio;
    g.x2o = -(((rx2 - x2 * sx) * nxo + (ry2 - x2 * sy) * nyo) * dsio);
    g.x2p = ((rx2 - x2 * sx) * nxp + (ry2 - x2 * sy) * nyp) * dsio - (sx * nxp + sy * nyp);
    g.yyo = ((rx1 + x1 * sy) * nyo - (ry1 - x1 * sx) * nxo) * dsio - (sx * nyo - sy * nxo);
    g.yyp = -(((rx1 - x1 * sy) * nyp - (ry1 + x1 * sx) * nxp) * dsio);
    return g;
}

void psilin(xfoil_t *xf, int i, double xi, double yi, double nxi, double nyi,
            double *psi, double *psi_ni, int geolin, int siglin)
{
    /* distance tolerance for determining if two points are the same */
    const double seps = (xf->s[xf->n - 1] - xf->s[0]) * 1.0e-5;

    const int io = i;

    xf->cosa = cos(xf->alfa);
    xf->sina = sin(xf->alfa);

    for (int jo = 0; jo < xf->n; jo++) {
        xf->dzdg[jo] = 0.0;
        xf->dzdn[jo] = 0.0;
        xf->dqdg[jo] = 0.0;
    }
    for (int jo = 0; jo < xf->n; jo++) {
        xf->dzdm[jo] = 0.0;
        xf->dqdm[jo] = 0.0;
    }

    xf->z_qinf = 0.0;
    xf->z_alfa = 0.0;
    xf->z_qdof0 = 0.0;
    xf->z_qdof1 = 0.0;
    xf->z_qdof2 = 0.0;
    xf->z_qdof3 = 0.0;

    *psi = 0.0;
    *psi_ni = 0.0;

    xf->qtan1 = 0.0;
    xf->qtan2 = 0.0;

    double scs, sds;
    if (xf->sharp) {
        scs = 1.0;
        sds = 0.0;
    } else {
        scs = xf->ante / xf->dste;
        sds = xf->aste / xf->dste;
    }

    /* Panel quantities: after the loop they still describe the TE panel
     * (the loop leaves with jo = n-1, jp = 0), which the TE terms use. */
    int jo, jp = 0;
    double apan = 0.0, x1 = 0.0, x2 = 0.0, yy = 0.0, g1 = 0.0, g2 = 0.0, t1 = 0.0, t2 = 0.0;
    double x1i = 0.0, x2i = 0.0, yyi = 0.0;
    struct geo gd = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    int closed_te = 0;
    for (jo = 0; jo < xf->n; jo++) {
        jp = jo + 1;
        int jm = jo - 1;
        int jq = jp + 1;

        if (jo == 0) {
            jm = jo;
        } else if (jo == xf->n - 2) {
            jq = jp;
        } else if (jo == xf->n - 1) {
            jp = 0;
            /* closed TE: no TE panel, go straight to the freestream terms */
            if (SQ(xf->x[jo] - xf->x[jp]) + SQ(xf->y[jo] - xf->y[jp]) < SQ(seps)) {
                closed_te = 1;
                break;
            }
        }

        const double dso = sqrt(SQ(xf->x[jo] - xf->x[jp]) + SQ(xf->y[jo] - xf->y[jp]));

        /* skip null panel */
        if (dso != 0.0) {
            const double dsio = 1.0 / dso;

            apan = xf->apanel[jo];

            const double rx1 = xi - xf->x[jo];
            const double ry1 = yi - xf->y[jo];
            const double rx2 = xi - xf->x[jp];
            const double ry2 = yi - xf->y[jp];

            const double sx = (xf->x[jp] - xf->x[jo]) * dsio;
            const double sy = (xf->y[jp] - xf->y[jo]) * dsio;

            x1 = sx * rx1 + sy * ry1;
            x2 = sx * rx2 + sy * ry2;
            yy = sx * ry1 - sy * rx1;

            const double rs1 = rx1 * rx1 + ry1 * ry1;
            const double rs2 = rx2 * rx2 + ry2 * ry2;

            /* reflection flag sgn avoids branch problems with arctan */
            const double sgn = (io >= 0 && io < xf->n) ? 1.0 : copysign(1.0, yy);

            /* log(r^2) and arctan(x/y), corrected for reflection */
            if (io != jo && rs1 > 0.0) {
                g1 = log(rs1);
                t1 = atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * xf->pi;
            } else {
                g1 = 0.0;
                t1 = 0.0;
            }
            if (io != jp && rs2 > 0.0) {
                g2 = log(rs2);
                t2 = atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * xf->pi;
            } else {
                g2 = 0.0;
                t2 = 0.0;
            }

            x1i = sx * nxi + sy * nyi;
            x2i = sx * nxi + sy * nyi;
            yyi = sx * nyi - sy * nxi;

            if (geolin) gd = geo_derivs(xf, jo, jp, rx1, ry1, rx2, ry2, x1, x2, sx, sy, dsio);

            if (jo == xf->n - 1) break;

            if (siglin)
                source_half_panels(xf, jm, jo, jp, jq, x1, x2, yy, rs1, rs2, g1, g2, t1, t2, apan, sgn, +1,
                                   x1i, x2i, yyi, dsio, psi, psi_ni);

            /* vortex panel contribution to Psi */
            const double dxinv = 1.0 / (x1 - x2);
            const double psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
            const double psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

            const double psx1 = 0.5 * g1;
            const double psx2 = -(0.5 * g2);
            const double psyy = t1 - t2;

            const double pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
            const double pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
            const double pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;

            const double gsum1 = xf->gamu[0][jp] + xf->gamu[0][jo];
            const double gsum2 = xf->gamu[1][jp] + xf->gamu[1][jo];
            const double gdif1 = xf->gamu[0][jp] - xf->gamu[0][jo];
            const double gdif2 = xf->gamu[1][jp] - xf->gamu[1][jo];

            const double gsum = xf->gam[jp] + xf->gam[jo];
            const double gdif = xf->gam[jp] - xf->gam[jo];

            *psi = *psi + xf->qopi * (psis * gsum + psid * gdif);

            /* dPsi/dGam */
            xf->dzdg[jo] = xf->dzdg[jo] + xf->qopi * (psis - psid);
            xf->dzdg[jp] = xf->dzdg[jp] + xf->qopi * (psis + psid);

            /* dPsi/dni */
            const double psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
            const double pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
            *psi_ni = *psi_ni + xf->qopi * (gsum * psni + gdif * pdni);

            xf->qtan1 = xf->qtan1 + xf->qopi * (gsum1 * psni + gdif1 * pdni);
            xf->qtan2 = xf->qtan2 + xf->qopi * (gsum2 * psni + gdif2 * pdni);

            xf->dqdg[jo] = xf->dqdg[jo] + xf->qopi * (psni - pdni);
            xf->dqdg[jp] = xf->dqdg[jp] + xf->qopi * (psni + pdni);

            if (geolin) {
                /* dPsi/dn */
                xf->dzdn[jo] = xf->dzdn[jo] + xf->qopi * gsum * (psx1 * gd.x1o + psx2 * gd.x2o + psyy * gd.yyo)
                           + xf->qopi * gdif * (pdx1 * gd.x1o + pdx2 * gd.x2o + pdyy * gd.yyo);
                xf->dzdn[jp] = xf->dzdn[jp] + xf->qopi * gsum * (psx1 * gd.x1p + psx2 * gd.x2p + psyy * gd.yyp)
                           + xf->qopi * gdif * (pdx1 * gd.x1p + pdx2 * gd.x2p + pdyy * gd.yyp);
                /* dPsi/dP */
                xf->z_qdof0 = xf->z_qdof0 + xf->qopi * ((psis - psid) * xf->qf0[jo] + (psis + psid) * xf->qf0[jp]);
                xf->z_qdof1 = xf->z_qdof1 + xf->qopi * ((psis - psid) * xf->qf1[jo] + (psis + psid) * xf->qf1[jp]);
                xf->z_qdof2 = xf->z_qdof2 + xf->qopi * ((psis - psid) * xf->qf2[jo] + (psis + psid) * xf->qf2[jp]);
                xf->z_qdof3 = xf->z_qdof3 + xf->qopi * ((psis - psid) * xf->qf3[jo] + (psis + psid) * xf->qf3[jp]);
            }
        }
    }

    if (!closed_te) {
        /* TE panel (jo = n-1, jp = 0) */
        const double psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
        const double pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

        const double psigx1 = -(t1 - apan);
        const double psigx2 = t2 - apan;
        const double psigyy = 0.5 * (g1 - g2);
        const double pgamx1 = 0.5 * g1;
        const double pgamx2 = -(0.5 * g2);
        const double pgamyy = t1 - t2;

        const double psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
        const double pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

        /* TE panel source and vortex strengths */
        const double sigte1 = 0.5 * scs * (xf->gamu[0][jp] - xf->gamu[0][jo]);
        const double sigte2 = 0.5 * scs * (xf->gamu[1][jp] - xf->gamu[1][jo]);
        const double gamte1 = -(0.5 * sds * (xf->gamu[0][jp] - xf->gamu[0][jo]));
        const double gamte2 = -(0.5 * sds * (xf->gamu[1][jp] - xf->gamu[1][jo]));

        xf->sigte = 0.5 * scs * (xf->gam[jp] - xf->gam[jo]);
        xf->gamte = -(0.5 * sds * (xf->gam[jp] - xf->gam[jo]));

        /* TE panel contribution to Psi */
        *psi = *psi + xf->hopi * (psig * xf->sigte + pgam * xf->gamte);

        /* dPsi/dGam */
        xf->dzdg[jo] = xf->dzdg[jo] - xf->hopi * psig * scs * 0.5;
        xf->dzdg[jp] = xf->dzdg[jp] + xf->hopi * psig * scs * 0.5;

        xf->dzdg[jo] = xf->dzdg[jo] + xf->hopi * pgam * sds * 0.5;
        xf->dzdg[jp] = xf->dzdg[jp] - xf->hopi * pgam * sds * 0.5;

        /* dPsi/dni */
        *psi_ni = *psi_ni + xf->hopi * (psigni * xf->sigte + pgamni * xf->gamte);

        xf->qtan1 = xf->qtan1 + xf->hopi * (psigni * sigte1 + pgamni * gamte1);
        xf->qtan2 = xf->qtan2 + xf->hopi * (psigni * sigte2 + pgamni * gamte2);

        xf->dqdg[jo] = xf->dqdg[jo] - xf->hopi * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
        xf->dqdg[jp] = xf->dqdg[jp] + xf->hopi * (psigni * 0.5 * scs - pgamni * 0.5 * sds);

        if (geolin) {
            /* dPsi/dn */
            xf->dzdn[jo] = xf->dzdn[jo] + xf->hopi * (psigx1 * gd.x1o + psigx2 * gd.x2o + psigyy * gd.yyo) * xf->sigte
                       + xf->hopi * (pgamx1 * gd.x1o + pgamx2 * gd.x2o + pgamyy * gd.yyo) * xf->gamte;
            xf->dzdn[jp] = xf->dzdn[jp] + xf->hopi * (psigx1 * gd.x1p + psigx2 * gd.x2p + psigyy * gd.yyp) * xf->sigte
                       + xf->hopi * (pgamx1 * gd.x1p + pgamx2 * gd.x2p + pgamyy * gd.yyp) * xf->gamte;

            /* dPsi/dP */
            xf->z_qdof0 = xf->z_qdof0 + xf->hopi * psig * 0.5 * (xf->qf0[jp] - xf->qf0[jo]) * scs
                      - xf->hopi * pgam * 0.5 * (xf->qf0[jp] - xf->qf0[jo]) * sds;
            xf->z_qdof1 = xf->z_qdof1 + xf->hopi * psig * 0.5 * (xf->qf1[jp] - xf->qf1[jo]) * scs
                      - xf->hopi * pgam * 0.5 * (xf->qf1[jp] - xf->qf1[jo]) * sds;
            xf->z_qdof2 = xf->z_qdof2 + xf->hopi * psig * 0.5 * (xf->qf2[jp] - xf->qf2[jo]) * scs
                      - xf->hopi * pgam * 0.5 * (xf->qf2[jp] - xf->qf2[jo]) * sds;
            xf->z_qdof3 = xf->z_qdof3 + xf->hopi * psig * 0.5 * (xf->qf3[jp] - xf->qf3[jo]) * scs
                      - xf->hopi * pgam * 0.5 * (xf->qf3[jp] - xf->qf3[jo]) * sds;
        }
    }

    /* freestream terms */
    *psi = *psi + xf->qinf * (xf->cosa * yi - xf->sina * xi);

    /* dPsi/dn */
    *psi_ni = *psi_ni + xf->qinf * (xf->cosa * nyi - xf->sina * nxi);

    xf->qtan1 = xf->qtan1 + xf->qinf * nyi;
    xf->qtan2 = xf->qtan2 - xf->qinf * nxi;

    /* dPsi/dQinf */
    xf->z_qinf = xf->z_qinf + (xf->cosa * yi - xf->sina * xi);

    /* dPsi/dalfa */
    xf->z_alfa = xf->z_alfa - xf->qinf * (xf->sina * yi + xf->cosa * xi);

    if (!xf->limage) return;

    /* ---- image system (ground effect) ---- */
    for (jo = 0; jo < xf->n; jo++) {
        jp = jo + 1;
        int jm = jo - 1;
        int jq = jp + 1;

        if (jo == 0) {
            jm = jo;
        } else if (jo == xf->n - 2) {
            jq = jp;
        } else if (jo == xf->n - 1) {
            jp = 0;
            if (SQ(xf->x[jo] - xf->x[jp]) + SQ(xf->y[jo] - xf->y[jp]) < SQ(seps)) return;
        }

        const double dso = sqrt(SQ(xf->x[jo] - xf->x[jp]) + SQ(xf->y[jo] - xf->y[jp]));

        /* skip null panel */
        if (dso != 0.0) {
            const double dsio = 1.0 / dso;

            apan = xf->pi - xf->apanel[jo] + 2.0 * xf->alfa;

            const double xjo = xf->x[jo] + 2.0 * (xf->yimage + xf->y[jo]) * xf->sina;
            const double yjo = xf->y[jo] - 2.0 * (xf->yimage + xf->y[jo]) * xf->cosa;
            const double xjp = xf->x[jp] + 2.0 * (xf->yimage + xf->y[jp]) * xf->sina;
            const double yjp = xf->y[jp] - 2.0 * (xf->yimage + xf->y[jp]) * xf->cosa;

            const double rx1 = xi - xjo;
            const double ry1 = yi - yjo;
            const double rx2 = xi - xjp;
            const double ry2 = yi - yjp;

            const double sx = (xjp - xjo) * dsio;
            const double sy = (yjp - yjo) * dsio;

            x1 = sx * rx1 + sy * ry1;
            x2 = sx * rx2 + sy * ry2;
            yy = sx * ry1 - sy * rx1;

            const double rs1 = rx1 * rx1 + ry1 * ry1;
            const double rs2 = rx2 * rx2 + ry2 * ry2;

            const double sgn = (io >= 0 && io < xf->n) ? 1.0 : copysign(1.0, yy);

            g1 = log(rs1);
            t1 = atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * xf->pi;

            g2 = log(rs2);
            t2 = atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * xf->pi;

            x1i = sx * nxi + sy * nyi;
            x2i = sx * nxi + sy * nyi;
            yyi = sx * nyi - sy * nxi;

            if (geolin) gd = geo_derivs(xf, jo, jp, rx1, ry1, rx2, ry2, x1, x2, sx, sy, dsio);

            if (jo == xf->n - 1) break;

            if (siglin)
                source_half_panels(xf, jm, jo, jp, jq, x1, x2, yy, rs1, rs2, g1, g2, t1, t2, apan, sgn, +1,
                                   x1i, x2i, yyi, dsio, psi, psi_ni);

            /* vortex panel contribution to Psi (image: opposite sign) */
            const double dxinv = 1.0 / (x1 - x2);
            const double psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
            const double psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

            const double psx1 = 0.5 * g1;
            const double psx2 = -(0.5 * g2);
            const double psyy = t1 - t2;

            const double pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
            const double pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
            const double pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;

            const double gsum1 = xf->gamu[0][jp] + xf->gamu[0][jo];
            const double gsum2 = xf->gamu[1][jp] + xf->gamu[1][jo];
            const double gdif1 = xf->gamu[0][jp] - xf->gamu[0][jo];
            const double gdif2 = xf->gamu[1][jp] - xf->gamu[1][jo];

            const double gsum = xf->gam[jp] + xf->gam[jo];
            const double gdif = xf->gam[jp] - xf->gam[jo];

            *psi = *psi - xf->qopi * (psis * gsum + psid * gdif);

            /* dPsi/dGam */
            xf->dzdg[jo] = xf->dzdg[jo] - xf->qopi * (psis - psid);
            xf->dzdg[jp] = xf->dzdg[jp] - xf->qopi * (psis + psid);

            /* dPsi/dni */
            const double psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
            const double pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
            *psi_ni = *psi_ni - xf->qopi * (gsum * psni + gdif * pdni);

            xf->qtan1 = xf->qtan1 - xf->qopi * (gsum1 * psni + gdif1 * pdni);
            xf->qtan2 = xf->qtan2 - xf->qopi * (gsum2 * psni + gdif2 * pdni);

            xf->dqdg[jo] = xf->dqdg[jo] - xf->qopi * (psni - pdni);
            xf->dqdg[jp] = xf->dqdg[jp] - xf->qopi * (psni + pdni);

            if (geolin) {
                /* dPsi/dn */
                xf->dzdn[jo] = xf->dzdn[jo] - xf->qopi * gsum * (psx1 * gd.x1o + psx2 * gd.x2o + psyy * gd.yyo)
                           - xf->qopi * gdif * (pdx1 * gd.x1o + pdx2 * gd.x2o + pdyy * gd.yyo);
                xf->dzdn[jp] = xf->dzdn[jp] - xf->qopi * gsum * (psx1 * gd.x1p + psx2 * gd.x2p + psyy * gd.yyp)
                           - xf->qopi * gdif * (pdx1 * gd.x1p + pdx2 * gd.x2p + pdyy * gd.yyp);
                /* dPsi/dP */
                xf->z_qdof0 = xf->z_qdof0 - xf->qopi * ((psis - psid) * xf->qf0[jo] + (psis + psid) * xf->qf0[jp]);
                xf->z_qdof1 = xf->z_qdof1 - xf->qopi * ((psis - psid) * xf->qf1[jo] + (psis + psid) * xf->qf1[jp]);
                xf->z_qdof2 = xf->z_qdof2 - xf->qopi * ((psis - psid) * xf->qf2[jo] + (psis + psid) * xf->qf2[jp]);
                xf->z_qdof3 = xf->z_qdof3 - xf->qopi * ((psis - psid) * xf->qf3[jo] + (psis + psid) * xf->qf3[jp]);
            }
        }
    }

    {
        /* image TE panel (jo = n-1, jp = 0) */
        const double psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
        const double pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

        const double psigx1 = -(t1 - apan);
        const double psigx2 = t2 - apan;
        const double psigyy = 0.5 * (g1 - g2);
        const double pgamx1 = 0.5 * g1;
        const double pgamx2 = -(0.5 * g2);
        const double pgamyy = t1 - t2;

        const double psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
        const double pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

        const double sigte1 = 0.5 * scs * (xf->gamu[0][jp] - xf->gamu[0][jo]);
        const double sigte2 = 0.5 * scs * (xf->gamu[1][jp] - xf->gamu[1][jo]);
        const double gamte1 = -(0.5 * sds * (xf->gamu[0][jp] - xf->gamu[0][jo]));
        const double gamte2 = -(0.5 * sds * (xf->gamu[1][jp] - xf->gamu[1][jo]));

        xf->sigte = 0.5 * scs * (xf->gam[jp] - xf->gam[jo]);
        xf->gamte = -(0.5 * sds * (xf->gam[jp] - xf->gam[jo]));

        *psi = *psi + xf->hopi * (psig * xf->sigte - pgam * xf->gamte);

        xf->dzdg[jo] = xf->dzdg[jo] - xf->hopi * psig * scs * 0.5;
        xf->dzdg[jp] = xf->dzdg[jp] + xf->hopi * psig * scs * 0.5;

        xf->dzdg[jo] = xf->dzdg[jo] - xf->hopi * pgam * sds * 0.5;
        xf->dzdg[jp] = xf->dzdg[jp] + xf->hopi * pgam * sds * 0.5;

        *psi_ni = *psi_ni + xf->hopi * (psigni * xf->sigte - pgamni * xf->gamte);

        xf->qtan1 = xf->qtan1 + xf->hopi * (psigni * sigte1 - pgamni * gamte1);
        xf->qtan2 = xf->qtan2 + xf->hopi * (psigni * sigte2 - pgamni * gamte2);

        xf->dqdg[jo] = xf->dqdg[jo] - xf->hopi * (psigni * 0.5 * scs + pgamni * 0.5 * sds);
        xf->dqdg[jp] = xf->dqdg[jp] + xf->hopi * (psigni * 0.5 * scs + pgamni * 0.5 * sds);

        if (geolin) {
            xf->dzdn[jo] = xf->dzdn[jo] + xf->hopi * (psigx1 * gd.x1o + psigx2 * gd.x2o + psigyy * gd.yyo) * xf->sigte
                       - xf->hopi * (pgamx1 * gd.x1o + pgamx2 * gd.x2o + pgamyy * gd.yyo) * xf->gamte;
            xf->dzdn[jp] = xf->dzdn[jp] + xf->hopi * (psigx1 * gd.x1p + psigx2 * gd.x2p + psigyy * gd.yyp) * xf->sigte
                       - xf->hopi * (pgamx1 * gd.x1p + pgamx2 * gd.x2p + pgamyy * gd.yyp) * xf->gamte;

            xf->z_qdof0 = xf->z_qdof0 + xf->hopi * psig * 0.5 * (xf->qf0[jp] - xf->qf0[jo]) * scs
                      + xf->hopi * pgam * 0.5 * (xf->qf0[jp] - xf->qf0[jo]) * sds;
            xf->z_qdof1 = xf->z_qdof1 + xf->hopi * psig * 0.5 * (xf->qf1[jp] - xf->qf1[jo]) * scs
                      + xf->hopi * pgam * 0.5 * (xf->qf1[jp] - xf->qf1[jo]) * sds;
            xf->z_qdof2 = xf->z_qdof2 + xf->hopi * psig * 0.5 * (xf->qf2[jp] - xf->qf2[jo]) * scs
                      + xf->hopi * pgam * 0.5 * (xf->qf2[jp] - xf->qf2[jo]) * sds;
            xf->z_qdof3 = xf->z_qdof3 + xf->hopi * psig * 0.5 * (xf->qf3[jp] - xf->qf3[jo]) * scs
                      + xf->hopi * pgam * 0.5 * (xf->qf3[jp] - xf->qf3[jo]) * sds;
        }
    }
}

void pswlin(xfoil_t *xf, int i, double xi, double yi, double nxi, double nyi, double *psi, double *psi_ni)
{
    const int io = i;

    xf->cosa = cos(xf->alfa);
    xf->sina = sin(xf->alfa);

    for (int jo = xf->n; jo < xf->n + xf->nw; jo++) {
        xf->dzdm[jo] = 0.0;
        xf->dqdm[jo] = 0.0;
    }

    *psi = 0.0;
    *psi_ni = 0.0;

    for (int jo = xf->n; jo < xf->n + xf->nw - 1; jo++) {
        const int jp = jo + 1;
        int jm = jo - 1;
        int jq = jp + 1;
        if (jo == xf->n)
            jm = jo;
        else if (jo == xf->n + xf->nw - 2)
            jq = jp;

        const double dso = sqrt(SQ(xf->x[jo] - xf->x[jp]) + SQ(xf->y[jo] - xf->y[jp]));
        const double dsio = 1.0 / dso;

        const double apan = xf->apanel[jo];

        const double rx1 = xi - xf->x[jo];
        const double ry1 = yi - xf->y[jo];
        const double rx2 = xi - xf->x[jp];
        const double ry2 = yi - xf->y[jp];

        const double sx = (xf->x[jp] - xf->x[jo]) * dsio;
        const double sy = (xf->y[jp] - xf->y[jo]) * dsio;

        const double x1 = sx * rx1 + sy * ry1;
        const double x2 = sx * rx2 + sy * ry2;
        const double yy = sx * ry1 - sy * rx1;

        const double rs1 = rx1 * rx1 + ry1 * ry1;
        const double rs2 = rx2 * rx2 + ry2 * ry2;

        /* on a wake node no reflection is needed */
        const double sgn = (io >= xf->n && io < xf->n + xf->nw) ? 1.0 : copysign(1.0, yy);

        double g1, t1, g2, t2;
        if (io != jo && rs1 > 0.0) {
            g1 = log(rs1);
            t1 = atan2(sgn * x1, sgn * yy) - (0.5 - 0.5 * sgn) * xf->pi;
        } else {
            g1 = 0.0;
            t1 = 0.0;
        }
        if (io != jp && rs2 > 0.0) {
            g2 = log(rs2);
            t2 = atan2(sgn * x2, sgn * yy) - (0.5 - 0.5 * sgn) * xf->pi;
        } else {
            g2 = 0.0;
            t2 = 0.0;
        }

        const double x1i = sx * nxi + sy * nyi;
        const double x2i = sx * nxi + sy * nyi;
        const double yyi = sx * nyi - sy * nxi;

        source_half_panels(xf, jm, jo, jp, jq, x1, x2, yy, rs1, rs2, g1, g2, t1, t2, apan, sgn, -1,
                           x1i, x2i, yyi, dsio, psi, psi_ni);
    }
}

/* ------------------------------------------------------------------------- */
/* Part 2: ggcalc, qwcalc, qdcalc. Fortran matrices keep their layout, so
 * AIJ(i,j) is aij[j][i], BIJ(i,j) bij[j][i], DIJ(i,j) dij[j][i], CIJ(iw,j)
 * cij[j][iw]; a Fortran column argument BIJ(1,j) is &bij[j][0]. */

void ggcalc(xfoil_t *xf)
{
    double psi, psi_n;

    /* distance of internal control point ahead of sharp TE
     * (fraction of smaller panel length adjacent to TE) */
    const double bwt = 0.1;

    if (xfoil_print) printf(" Calculating unit vorticity distributions ...\n");

    for (int i = 0; i < xf->n; i++) {
        xf->gam[i] = 0.0;
        xf->gamu[0][i] = 0.0;
        xf->gamu[1][i] = 0.0;
    }

    /* matrix system for Psi = Psio on the airfoil surface;
     * unknowns are (dGamma)i and dPsio */
    for (int i = 0; i < xf->n; i++) {
        /* Psi and dPsi/dGamma for the current node */
        psilin(xf, i, xf->x[i], xf->y[i], xf->nx[i], xf->ny[i], &psi, &psi_n, 0, 1);

        /* RES1 = PSI( 0) - PSIO, RES2 = PSI(90) - PSIO */
        const double res1 = xf->qinf * xf->y[i];
        const double res2 = -(xf->qinf * xf->x[i]);

        /* dRes/dGamma */
        for (int j = 0; j < xf->n; j++) xf->aij[j][i] = xf->dzdg[j];
        for (int j = 0; j < xf->n; j++) xf->bij[j][i] = -xf->dzdm[j];

        /* dRes/dPsio */
        xf->aij[xf->n][i] = -1.0;

        xf->gamu[0][i] = -res1;
        xf->gamu[1][i] = -res2;
    }

    /* Kutta condition: RES = GAM(1) + GAM(N) */
    const double res = 0.0;

    for (int j = 0; j < xf->n + 1; j++) xf->aij[j][xf->n] = 0.0;

    xf->aij[0][xf->n] = 1.0;
    xf->aij[xf->n - 1][xf->n] = 1.0;

    xf->gamu[0][xf->n] = -res;
    xf->gamu[1][xf->n] = -res;

    /* Kutta condition has no direct source influence */
    for (int j = 0; j < xf->n; j++) xf->bij[j][xf->n] = 0.0;

    if (xf->sharp) {
        /* zero internal velocity in the TE corner: TE bisector angle */
        const double ag1 = atan2(-xf->yp[0], -xf->xp[0]);
        const double ag2 = atanc(xf->yp[xf->n - 1], xf->xp[xf->n - 1], ag1);
        const double abis = 0.5 * (ag1 + ag2);
        const double cbis = cos(abis);
        const double sbis = sin(abis);

        /* minimum panel length adjacent to the TE */
        const double ds1 = sqrt(SQ(xf->x[0] - xf->x[1]) + SQ(xf->y[0] - xf->y[1]));
        const double ds2 = sqrt(SQ(xf->x[xf->n - 1] - xf->x[xf->n - 2]) + SQ(xf->y[xf->n - 1] - xf->y[xf->n - 2]));
        const double dsmin = fmin_f(ds1, ds2);

        /* control point on the bisector just ahead of the TE point */
        const double xbis = xf->xte - bwt * dsmin * cbis;
        const double ybis = xf->yte - bwt * dsmin * sbis;

        /* velocity component along the bisector line (not a panel node: -1) */
        double qbis;
        psilin(xf, -1, xbis, ybis, -sbis, cbis, &psi, &qbis, 0, 1);

        /* dRes/dGamma, -dRes/dMass, dRes/dPsio (RES = qbis) */
        for (int j = 0; j < xf->n; j++) xf->aij[j][xf->n - 1] = xf->dqdg[j];
        for (int j = 0; j < xf->n; j++) xf->bij[j][xf->n - 1] = -xf->dqdm[j];
        xf->aij[xf->n][xf->n - 1] = 0.0;

        /* -dRes/dUinf, -dRes/dVinf */
        xf->gamu[0][xf->n - 1] = -cbis;
        xf->gamu[1][xf->n - 1] = -sbis;
    }

    /* LU-factor the coefficient matrix AIJ */
    ludcmp(IQX, xf->n + 1, &xf->aij[0][0], xf->aijpiv);
    xf->lqaij = 1;

    /* solve for the two vorticity distributions */
    baksub(IQX, xf->n + 1, &xf->aij[0][0], xf->aijpiv, &xf->gamu[0][0]);
    baksub(IQX, xf->n + 1, &xf->aij[0][0], xf->aijpiv, &xf->gamu[1][0]);

    /* inviscid alpha = 0, 90 surface speeds for this geometry */
    for (int i = 0; i < xf->n; i++) {
        xf->qinvu[0][i] = xf->gamu[0][i];
        xf->qinvu[1][i] = xf->gamu[1][i];
    }

    xf->lgamu = 1;
}

void qwcalc(xfoil_t *xf)
{
    double psi, psi_ni;

    /* first wake point (same as TE) */
    xf->qinvu[0][xf->n] = xf->qinvu[0][xf->n - 1];
    xf->qinvu[1][xf->n] = xf->qinvu[1][xf->n - 1];

    /* rest of the wake */
    for (int i = xf->n + 1; i < xf->n + xf->nw; i++) {
        psilin(xf, i, xf->x[i], xf->y[i], xf->nx[i], xf->ny[i], &psi, &psi_ni, 0, 0);
        xf->qinvu[0][i] = xf->qtan1;
        xf->qinvu[1][i] = xf->qtan2;
    }
}

void qdcalc(xfoil_t *xf)
{
    double psi, psi_n;

    if (xfoil_print) printf(" Calculating source influence matrix ...\n");

    if (!xf->ladij) {
        /* source influence matrix for the airfoil surface, if it doesn't exist */
        for (int j = 0; j < xf->n; j++) {
            /* multiply each dPsi/Sig vector by the inverse of the factored dPsi/dGam matrix */
            baksub(IQX, xf->n + 1, &xf->aij[0][0], xf->aijpiv, &xf->bij[j][0]);

            /* store the resulting dGam/dSig = dQtan/dSig vector */
            for (int i = 0; i < xf->n; i++) xf->dij[j][i] = xf->bij[j][i];
        }
        xf->ladij = 1;
    }

    /* coefficient matrix of dPsi/dm on the airfoil surface */
    for (int i = 0; i < xf->n; i++) {
        pswlin(xf, i, xf->x[i], xf->y[i], xf->nx[i], xf->ny[i], &psi, &psi_n);
        for (int j = xf->n; j < xf->n + xf->nw; j++) xf->bij[j][i] = -xf->dzdm[j];
    }

    /* Kutta condition (no direct source influence) */
    for (int j = xf->n; j < xf->n + xf->nw; j++) xf->bij[j][xf->n] = 0.0;

    /* sharp TE gamma extrapolation also has no source influence */
    if (xf->sharp) {
        for (int j = xf->n; j < xf->n + xf->nw; j++) xf->bij[j][xf->n - 1] = 0.0;
    }

    /* multiply by the inverse of the factored dPsi/dGam matrix */
    for (int j = xf->n; j < xf->n + xf->nw; j++) baksub(IQX, xf->n + 1, &xf->aij[0][0], xf->aijpiv, &xf->bij[j][0]);

    /* source influence matrix for the wake sources */
    for (int i = 0; i < xf->n; i++)
        for (int j = xf->n; j < xf->n + xf->nw; j++) xf->dij[j][i] = xf->bij[j][i];

    /* influence of the sources on the wake velocities:
     * dQtan/dGam and dQtan/dSig at the wake points */
    for (int i = xf->n; i < xf->n + xf->nw; i++) {
        const int iw = i - xf->n;       /* 0-based wake node (Fortran IW = i - N, 1-based) */
        if (iw >= IWX || iw < 0) break;

        /* airfoil contribution at the wake panel node */
        psilin(xf, i, xf->x[i], xf->y[i], xf->nx[i], xf->ny[i], &psi, &psi_n, 0, 1);

        for (int j = 0; j < xf->n; j++) xf->cij[j][iw] = xf->dqdg[j];
        for (int j = 0; j < xf->n; j++) xf->dij[j][i] = xf->dqdm[j];

        /* wake contribution */
        pswlin(xf, i, xf->x[i], xf->y[i], xf->nx[i], xf->ny[i], &psi, &psi_n);

        for (int j = xf->n; j < xf->n + xf->nw; j++) xf->dij[j][i] = xf->dqdm[j];
    }

    /* add the effect of all sources on the airfoil vorticity, which affects wake Qtan */
    for (int i = xf->n; i < xf->n + xf->nw; i++) {
        const int iw = i - xf->n;

        /* airfoil surface source contribution first */
        for (int j = 0; j < xf->n; j++) {
            double sum = 0.0;
            for (int k = 0; k < xf->n; k++) sum = sum + xf->cij[k][iw] * xf->dij[j][k];
            xf->dij[j][i] = xf->dij[j][i] + sum;
        }

        /* wake source contribution next */
        for (int j = xf->n; j < xf->n + xf->nw; j++) {
            double sum = 0.0;
            for (int k = 0; k < xf->n; k++) sum = sum + xf->cij[k][iw] * xf->bij[j][k];
            xf->dij[j][i] = xf->dij[j][i] + sum;
        }
    }

    /* the first wake point has the same velocity as the trailing edge */
    for (int j = 0; j < xf->n + xf->nw; j++) xf->dij[j][xf->n] = xf->dij[j][xf->n - 1];

    xf->lwdij = 1;
}

/* ------------------------------------------------------------------------- */
/* Part 3: xywake, stfind, iblpan, xicalc. ist, iblte and ipan hold 0-based
 * indices (Fortran value - 1); nbl holds counts, as in the Fortran. */

void xywake(xfoil_t *xf)
{
    double psi, psi_x, psi_y;

    if (xfoil_print) printf(" Calculating wake trajectory ...\n");

    /* number of wake points */
    xf->nw = xf->n / 12 + 10 * (int)xf->waklen;
    if (xf->nw > IWX) {
        if (xfoil_print) printf(" Array size (IWX) too small.  Last wake point index reduced.\n");
        xf->nw = IWX;
    }

    const double ds1 = 0.5 * (xf->s[1] - xf->s[0] + xf->s[xf->n - 1] - xf->s[xf->n - 2]);
    setexp(&xf->snew[xf->n], ds1, xf->waklen * xf->chord, xf->nw);

    xf->xte = 0.5 * (xf->x[0] + xf->x[xf->n - 1]);
    xf->yte = 0.5 * (xf->y[0] + xf->y[xf->n - 1]);

    /* first wake point a tiny distance behind the TE */
    int i = xf->n;
    const double sx = 0.5 * (xf->yp[xf->n - 1] - xf->yp[0]);
    const double sy = 0.5 * (xf->xp[0] - xf->xp[xf->n - 1]);
    const double smod = sqrt(SQ(sx) + SQ(sy));
    xf->nx[i] = sx / smod;
    xf->ny[i] = sy / smod;
    xf->x[i] = xf->xte - 0.0001 * xf->ny[i];
    xf->y[i] = xf->yte + 0.0001 * xf->nx[i];
    xf->s[i] = xf->s[xf->n - 1];

    /* streamfunction gradient components at the first point */
    psilin(xf, i, xf->x[i], xf->y[i], 1.0, 0.0, &psi, &psi_x, 0, 0);
    psilin(xf, i, xf->x[i], xf->y[i], 0.0, 1.0, &psi, &psi_y, 0, 0);

    /* unit vector normal to the wake at the first point */
    xf->nx[i + 1] = -(psi_x / sqrt(SQ(psi_x) + SQ(psi_y)));
    xf->ny[i + 1] = -(psi_y / sqrt(SQ(psi_x) + SQ(psi_y)));

    /* angle of the wake panel normal */
    xf->apanel[i] = atan2(psi_y, psi_x);

    /* rest of the wake points */
    for (i = xf->n + 1; i < xf->n + xf->nw; i++) {
        const double ds = xf->snew[i] - xf->snew[i - 1];

        /* new point ds downstream of the last point */
        xf->x[i] = xf->x[i - 1] - ds * xf->ny[i];
        xf->y[i] = xf->y[i - 1] + ds * xf->nx[i];
        xf->s[i] = xf->s[i - 1] + ds;

        if (i != xf->n + xf->nw - 1) {
            /* normal vector for the next point */
            psilin(xf, i, xf->x[i], xf->y[i], 1.0, 0.0, &psi, &psi_x, 0, 0);
            psilin(xf, i, xf->x[i], xf->y[i], 0.0, 1.0, &psi, &psi_y, 0, 0);

            xf->nx[i + 1] = -(psi_x / sqrt(SQ(psi_x) + SQ(psi_y)));
            xf->ny[i + 1] = -(psi_y / sqrt(SQ(psi_x) + SQ(psi_y)));

            xf->apanel[i] = atan2(psi_y, psi_x);
        }
    }

    /* wake presence flag and corresponding alpha */
    xf->lwake = 1;
    xf->awake = xf->alfa;

    /* the old source influence matrix is invalid for the new wake geometry */
    xf->lwdij = 0;
}

void stfind(xfoil_t *xf)
{
    /* first panel with a sign change of gam (+ to -) */
    int i = 0;
    while (i < xf->n - 1 && !(xf->gam[i] >= 0.0 && xf->gam[i + 1] < 0.0)) i++;

    if (i >= xf->n - 1) {
        if (xfoil_print) printf(" STFIND: Stagnation point not found. Continuing ...\n");
        i = xf->n / 2 - 1;      /* Fortran I = N/2 (1-based) */
    }

    xf->ist = i;
    const double dgam = xf->gam[i + 1] - xf->gam[i];
    const double ds = xf->s[i + 1] - xf->s[i];

    /* evaluate so as to minimize roundoff for very small GAM(I) or GAM(I+1) */
    if (xf->gam[i] < -xf->gam[i + 1])
        xf->sst = xf->s[i] - ds * (xf->gam[i] / dgam);
    else
        xf->sst = xf->s[i + 1] - ds * (xf->gam[i + 1] / dgam);

    /* tweak the stagnation point if it falls right on a node (very unlikely) */
    if (xf->sst <= xf->s[i]) xf->sst = xf->s[i] + 1.0e-7;
    if (xf->sst >= xf->s[i + 1]) xf->sst = xf->s[i + 1] - 1.0e-7;

    xf->sst_go = (xf->sst - xf->s[i + 1]) / dgam;
    xf->sst_gp = (xf->s[i] - xf->sst) / dgam;
}

void iblpan(xfoil_t *xf)
{
    int ibl;

    /* top surface first: stations 1.. go from ist back to node 0 */
    ibl = 0;
    for (int i = xf->ist; i >= 0; i--) {
        ibl = ibl + 1;
        xf->ipan[0][ibl] = i;
        xf->vti[0][ibl] = 1.0;
    }
    xf->iblte[0] = ibl;
    xf->nbl[0] = ibl + 1;                   /* count (Fortran NBL = last 1-based station) */

    /* bottom surface next */
    ibl = 0;
    for (int i = xf->ist + 1; i < xf->n; i++) {
        ibl = ibl + 1;
        xf->ipan[1][ibl] = i;
        xf->vti[1][ibl] = -1.0;
    }

    /* wake */
    xf->iblte[1] = ibl;

    for (int iw = 1; iw <= xf->nw; iw++) {
        const int i = xf->n + iw - 1;       /* Fortran I = N + IW */
        ibl = xf->iblte[1] + iw;
        xf->ipan[1][ibl] = i;
        xf->vti[1][ibl] = -1.0;
    }

    xf->nbl[1] = xf->iblte[1] + xf->nw + 1;

    /* upper wake pointers (for plotting only) */
    for (int iw = 1; iw <= xf->nw; iw++) {
        xf->ipan[0][xf->iblte[0] + iw] = xf->ipan[1][xf->iblte[1] + iw];
        xf->vti[0][xf->iblte[0] + iw] = 1.0;
    }

    /* Fortran IBLMAX = max(IBLTE(1), IBLTE(2)) + NW with 1-based IBLTE */
    const int iblmax = (xf->iblte[0] > xf->iblte[1] ? xf->iblte[0] : xf->iblte[1]) + 1 + xf->nw;
    if (iblmax > IVX) {
        if (xfoil_print) {
            /* list-directed: a leading blank, then the literal (which starts with a blank), then I12 */
            printf("  ***  BL array overflow.\n");
            printf("  ***  Increase IVX to at least%12d\n", iblmax);
        }
        xfoil_stop("IBLPAN: BL array overflow (increase IVX)");   /* Fortran: bare stop */
    }

    xf->lipan = 1;
}

void xicalc(xfoil_t *xf)
{
    static const double xfeps = 1.0e-7;

    /* minimum xi node arc length near the stagnation point */
    const double xeps = xfeps * (xf->s[xf->n - 1] - xf->s[0]);

    xf->xssi[0][0] = 0.0;
    for (int ibl = 1; ibl <= xf->iblte[0]; ibl++) {
        const int i = xf->ipan[0][ibl];
        xf->xssi[0][ibl] = fmax_f(xf->sst - xf->s[i], xeps);
    }

    xf->xssi[1][0] = 0.0;
    for (int ibl = 1; ibl <= xf->iblte[1]; ibl++) {
        const int i = xf->ipan[1][ibl];
        xf->xssi[1][ibl] = fmax_f(xf->s[i] - xf->sst, xeps);
    }

    int ibl1 = xf->iblte[0] + 1;
    xf->xssi[0][ibl1] = xf->xssi[0][ibl1 - 1];

    int ibl2 = xf->iblte[1] + 1;
    xf->xssi[1][ibl2] = xf->xssi[1][ibl2 - 1];

    /* wake (side 2 stations; mirrored onto side 1) */
    for (int ibl = xf->iblte[1] + 2; ibl <= xf->nbl[1] - 1; ibl++) {
        const int i = xf->ipan[1][ibl];
        const double dxssi = sqrt(SQ(xf->x[i] - xf->x[i - 1]) + SQ(xf->y[i] - xf->y[i - 1]));

        ibl1 = xf->iblte[0] + ibl - xf->iblte[1];
        ibl2 = xf->iblte[1] + ibl - xf->iblte[1];
        xf->xssi[0][ibl1] = xf->xssi[0][ibl1 - 1] + dxssi;
        xf->xssi[1][ibl2] = xf->xssi[1][ibl2 - 1] + dxssi;
    }

    /* trailing edge flap length to TE gap ratio */
    const double telrat = 2.50;

    /* parameters for the TE flap cubics */
    const double crosp = (xf->xp[0] * xf->yp[xf->n - 1] - xf->yp[0] * xf->xp[xf->n - 1])
                         / sqrt((SQ(xf->xp[0]) + SQ(xf->yp[0])) * (SQ(xf->xp[xf->n - 1]) + SQ(xf->yp[xf->n - 1])));
    double dwdxte = crosp / sqrt(1.0 - SQ(crosp));

    /* limit the cubic to avoid absurd TE gap widths */
    dwdxte = fmax_f(dwdxte, -(3.0 / telrat));
    dwdxte = fmin_f(dwdxte, 3.0 / telrat);

    const double aa = 3.0 + telrat * dwdxte;
    const double bb = -2.0 - telrat * dwdxte;

    if (xf->sharp) {
        for (int iw = 0; iw < xf->nw; iw++) xf->wgap[iw] = 0.0;
    } else {
        /* TE flap (wake gap) array */
        for (int iw = 1; iw <= xf->nw; iw++) {
            const int ibl = xf->iblte[1] + iw;
            const double zn = 1.0 - (xf->xssi[1][ibl] - xf->xssi[1][xf->iblte[1]]) / (telrat * xf->ante);
            xf->wgap[iw - 1] = 0.0;
            if (zn >= 0.0) xf->wgap[iw - 1] = xf->ante * (aa + bb * zn) * SQ(zn);
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Part 4: uicalc, qvfue, qiset, gamqv, stmove, ueset. BL arrays such as
 * UEDG(ibl,is) are uedg[is][ibl]; Fortran station loops `do ibl = 2, NBL(is)`
 * become `for (ibl = 1; ibl < nbl[is]; ibl++)`. */

void uicalc(xfoil_t *xf)
{
    for (int is = 0; is < 2; is++) {
        xf->uinv[is][0] = 0.0;
        xf->uinv_a[is][0] = 0.0;
        for (int ibl = 1; ibl < xf->nbl[is]; ibl++) {
            const int i = xf->ipan[is][ibl];
            xf->uinv[is][ibl] = xf->vti[is][ibl] * xf->qinv[i];
            xf->uinv_a[is][ibl] = xf->vti[is][ibl] * xf->qinv_a[i];
        }
    }
}

void qvfue(xfoil_t *xf)
{
    for (int is = 0; is < 2; is++) {
        for (int ibl = 1; ibl < xf->nbl[is]; ibl++) {
            const int i = xf->ipan[is][ibl];
            xf->qvis[i] = xf->vti[is][ibl] * xf->uedg[is][ibl];
        }
    }
}

void qiset(xfoil_t *xf)
{
    xf->cosa = cos(xf->alfa);
    xf->sina = sin(xf->alfa);

    for (int i = 0; i < xf->n + xf->nw; i++) {
        xf->qinv[i] = xf->cosa * xf->qinvu[0][i] + xf->sina * xf->qinvu[1][i];
        xf->qinv_a[i] = -(xf->sina * xf->qinvu[0][i]) + xf->cosa * xf->qinvu[1][i];
    }
}

void gamqv(xfoil_t *xf)
{
    for (int i = 0; i < xf->n; i++) {
        xf->gam[i] = xf->qvis[i];
        xf->gam_a[i] = xf->qinv_a[i];
    }
}

void stmove(xfoil_t *xf)
{
    /* locate the new stagnation point arc length sst from the gam distribution */
    const int istold = xf->ist;
    stfind(xf);

    if (istold == xf->ist) {
        /* recalculate the arc length array */
        xicalc(xf);
    } else {
        /* new BL position -> panel position pointers, inviscid Ue, arc lengths,
         * BL position -> system line pointers */
        iblpan(xf);
        uicalc(xf);
        xicalc(xf);
        iblsys(xf);

        /* In the loops below b is the 0-based station (Fortran ibl - 1):
         * Fortran `do ibl = NBL, idif+2, -1` is b = nbl-1 .. idif+1,
         * `do ibl = idif+1, 2, -1` is b = idif .. 1, `do ibl = 2, NBL` is b = 1 .. nbl-1. */
        if (xf->ist > istold) {
            /* more points on the top side (is = 1) */
            const int idif = xf->ist - istold;

            xf->itran[0] = xf->itran[0] + idif;
            xf->itran[1] = xf->itran[1] - idif;

            /* move top side BL variables downstream */
            for (int b = xf->nbl[0] - 1; b >= idif + 1; b--) {
                xf->ctau[0][b] = xf->ctau[0][b - idif];
                xf->thet[0][b] = xf->thet[0][b - idif];
                xf->dstr[0][b] = xf->dstr[0][b - idif];
                xf->uedg[0][b] = xf->uedg[0][b - idif];
            }

            /* BL variables between the old and new stagnation point */
            const double dudx = xf->uedg[0][idif + 1] / xf->xssi[0][idif + 1];
            for (int b = idif; b >= 1; b--) {
                xf->ctau[0][b] = xf->ctau[0][idif + 1];
                xf->thet[0][b] = xf->thet[0][idif + 1];
                xf->dstr[0][b] = xf->dstr[0][idif + 1];
                xf->uedg[0][b] = dudx * xf->xssi[0][b];
            }

            /* move bottom side BL variables upstream */
            for (int b = 1; b <= xf->nbl[1] - 1; b++) {
                xf->ctau[1][b] = xf->ctau[1][b + idif];
                xf->thet[1][b] = xf->thet[1][b + idif];
                xf->dstr[1][b] = xf->dstr[1][b + idif];
                xf->uedg[1][b] = xf->uedg[1][b + idif];
            }
        } else {
            /* more points on the bottom side (is = 2) */
            const int idif = istold - xf->ist;

            xf->itran[0] = xf->itran[0] - idif;
            xf->itran[1] = xf->itran[1] + idif;

            /* move bottom side BL variables downstream */
            for (int b = xf->nbl[1] - 1; b >= idif + 1; b--) {
                xf->ctau[1][b] = xf->ctau[1][b - idif];
                xf->thet[1][b] = xf->thet[1][b - idif];
                xf->dstr[1][b] = xf->dstr[1][b - idif];
                xf->uedg[1][b] = xf->uedg[1][b - idif];
            }

            /* BL variables between the old and new stagnation point */
            const double dudx = xf->uedg[1][idif + 1] / xf->xssi[1][idif + 1];
            for (int b = idif; b >= 1; b--) {
                xf->ctau[1][b] = xf->ctau[1][idif + 1];
                xf->thet[1][b] = xf->thet[1][idif + 1];
                xf->dstr[1][b] = xf->dstr[1][idif + 1];
                xf->uedg[1][b] = dudx * xf->xssi[1][b];
            }

            /* move top side BL variables upstream */
            for (int b = 1; b <= xf->nbl[0] - 1; b++) {
                xf->ctau[0][b] = xf->ctau[0][b + idif];
                xf->thet[0][b] = xf->thet[0][b + idif];
                xf->dstr[0][b] = xf->dstr[0][b + idif];
                xf->uedg[0][b] = xf->uedg[0][b + idif];
            }
        }

        /* tweak Ue so it's not zero, in case the stagnation point is right on a node */
        const double ueps = 1.0e-7;
        for (int is = 0; is < 2; is++) {
            for (int ibl = 1; ibl < xf->nbl[is]; ibl++) {
                const int i = xf->ipan[is][ibl];
                if (xf->uedg[is][ibl] <= ueps) {
                    xf->uedg[is][ibl] = ueps;
                    xf->qvis[i] = xf->vti[is][ibl] * ueps;
                    xf->gam[i] = xf->vti[is][ibl] * ueps;
                }
            }
        }
    }

    /* new mass array since Ue has been tweaked */
    for (int is = 0; is < 2; is++)
        for (int ibl = 1; ibl < xf->nbl[is]; ibl++) xf->mass[is][ibl] = xf->dstr[is][ibl] * xf->uedg[is][ibl];
}

void ueset(xfoil_t *xf)
{
    for (int is = 0; is < 2; is++) {
        for (int ibl = 1; ibl < xf->nbl[is]; ibl++) {
            const int i = xf->ipan[is][ibl];

            double dui = 0.0;
            for (int js = 0; js < 2; js++) {
                for (int jbl = 1; jbl < xf->nbl[js]; jbl++) {
                    const int j = xf->ipan[js][jbl];
                    const double ue_m = -(xf->vti[is][ibl] * xf->vti[js][jbl] * xf->dij[j][i]);
                    dui = dui + ue_m * xf->mass[js][jbl];
                }
            }

            xf->uedg[is][ibl] = xf->uinv[is][ibl] + dui;
        }
    }
}
