/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of m_xblsys.f90. Part 1: closure correlations (hkin, dil, dilw, hsl,
 *  cfl, hst, cft, hct). Fortran unary minus binds weaker than * and /, so
 *  `-a*b` is written -(a*b); x**n (integer) is __builtin_powi, x**r is pow().
 */
#include "xblsys.h"

#include <math.h>
#include <stdio.h>

#include "xstate.h"
#include "xfortran.h"

#define P2(v) __builtin_powi((v), 2)
#define P3(v) __builtin_powi((v), 3)

void hkin(double h, double msq, double *hk, double *hk_h, double *hk_msq)
{
    /* kinematic shape parameter (assuming air), from Whitfield */
    *hk = (h - 0.29 * msq) / (1.0 + 0.113 * msq);
    *hk_h = 1.0 / (1.0 + 0.113 * msq);
    *hk_msq = (-0.29 - 0.113 * *hk) / (1.0 + 0.113 * msq);
}

void dil(double hk, double rt, double *di, double *di_hk, double *di_rt)
{
    /* laminar dissipation function (2 CD/H*), from Falkner-Skan */
    if (hk < 4.0) {
        *di = (0.00205 * pow(4.0 - hk, 5.5) + 0.207) / rt;
        *di_hk = (-(0.00205 * 5.5 * pow(4.0 - hk, 4.5))) / rt;
    } else {
        const double hkb = hk - 4.0;
        const double den = 1.0 + 0.02 * P2(hkb);
        *di = (-(0.0016 * P2(hkb) / den) + 0.207) / rt;
        *di_hk = (-(0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * P2(hkb) / P2(den)))) / rt;
    }
    *di_rt = -(*di / rt);
}

void dilw(double hk, double rt, double *di, double *di_hk, double *di_rt)
{
    double hs, hs_hk, hs_rt, hs_msq;
    const double msq = 0.0;
    hsl(hk, rt, msq, &hs, &hs_hk, &hs_rt, &hs_msq);

    /* laminar wake dissipation function (2 CD/H*) */
    const double rcd = 1.10 * P2(1.0 - 1.0 / hk) / hk;
    const double rcd_hk = -(1.10 * (1.0 - 1.0 / hk) * 2.0 / P3(hk)) - rcd / hk;

    *di = 2.0 * rcd / (hs * rt);
    *di_hk = 2.0 * rcd_hk / (hs * rt) - (*di / hs) * hs_hk;
    *di_rt = -(*di / rt) - (*di / hs) * hs_rt;
}

void hsl(double hk, double rt, double msq, double *hs, double *hs_hk, double *hs_rt, double *hs_msq)
{
    (void)rt;
    (void)msq;
    /* laminar H* correlation */
    if (hk < 4.35) {
        const double tmp = hk - 4.35;
        *hs = 0.0111 * P2(tmp) / (hk + 1.0) - 0.0278 * P3(tmp) / (hk + 1.0) + 1.528 - 0.0002 * P2(tmp * hk);
        *hs_hk = 0.0111 * (2.0 * tmp - P2(tmp) / (hk + 1.0)) / (hk + 1.0)
                 - 0.0278 * (3.0 * P2(tmp) - P3(tmp) / (hk + 1.0)) / (hk + 1.0)
                 - 0.0002 * 2.0 * tmp * hk * (tmp + hk);
    } else {
        const double hs2 = 0.015;
        *hs = hs2 * P2(hk - 4.35) / hk + 1.528;
        *hs_hk = hs2 * 2.0 * (hk - 4.35) / hk - hs2 * P2(hk - 4.35) / P2(hk);
    }

    *hs_rt = 0.0;
    *hs_msq = 0.0;
}

void cfl(double hk, double rt, double msq, double *cf, double *cf_hk, double *cf_rt, double *cf_msq)
{
    (void)msq;
    /* laminar skin friction function (Cf), from Falkner-Skan */
    if (hk < 5.5) {
        const double tmp = P3(5.5 - hk) / (hk + 1.0);
        *cf = (0.0727 * tmp - 0.07) / rt;
        *cf_hk = (-(0.0727 * tmp * 3.0 / (5.5 - hk)) - 0.0727 * tmp / (hk + 1.0)) / rt;
    } else {
        const double tmp = 1.0 - 1.0 / (hk - 4.5);
        *cf = (0.015 * P2(tmp) - 0.07) / rt;
        *cf_hk = (0.015 * tmp * 2.0 / P2(hk - 4.5)) / rt;
    }
    *cf_rt = -(*cf / rt);
    *cf_msq = 0.0;
}

void hst(double hk, double rt, double msq, double *hs, double *hs_hk, double *hs_rt, double *hs_msq)
{
    static const double hsmin = 1.500, dhsinf = 0.015;

    /* turbulent H* correlation; limited Rtheta dependence for Rtheta < 200 (12/4/94) */
    double ho, ho_rt;
    if (rt > 400.0) {
        ho = 3.0 + 400.0 / rt;
        ho_rt = -(400.0 / P2(rt));
    } else {
        ho = 4.0;
        ho_rt = 0.0;
    }

    double rtz, rtz_rt;
    if (rt > 200.0) {
        rtz = rt;
        rtz_rt = 1.0;
    } else {
        rtz = 200.0;
        rtz_rt = 0.0;
    }

    if (hk < ho) {
        /* attached branch: new correlation 29 Nov 91 (arctan(y+) + Schlichting profiles) */
        const double hr = (ho - hk) / (ho - 1.0);
        const double hr_hk = -(1.0 / (ho - 1.0));
        const double hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;
        *hs = (2.0 - hsmin - 4.0 / rtz) * P2(hr) * 1.5 / (hk + 0.5) + hsmin + 4.0 / rtz;
        *hs_hk = -((2.0 - hsmin - 4.0 / rtz) * P2(hr) * 1.5 / P2(hk + 0.5))
                 + (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk;
        *hs_rt = (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt
                 + (P2(hr) * 1.5 / (hk + 0.5) - 1.0) * 4.0 / P2(rtz) * rtz_rt;
    } else {
        /* separated branch */
        const double grt = log(rtz);
        const double hdif = hk - ho;
        const double rtmp = hk - ho + 4.0 / grt;
        const double htmp = 0.007 * grt / P2(rtmp) + dhsinf / hk;
        const double htmp_hk = -(0.014 * grt / P3(rtmp)) - dhsinf / P2(hk);
        const double htmp_rt = -(0.014 * grt / P3(rtmp) * (-ho_rt - 4.0 / P2(grt) / rtz * rtz_rt))
                               + 0.007 / P2(rtmp) / rtz * rtz_rt;
        *hs = P2(hdif) * htmp + hsmin + 4.0 / rtz;
        *hs_hk = hdif * 2.0 * htmp + P2(hdif) * htmp_hk;
        *hs_rt = P2(hdif) * htmp_rt - 4.0 / P2(rtz) * rtz_rt + hdif * 2.0 * htmp * (-ho_rt);
    }

    /* Whitfield's minor additional compressibility correction */
    const double fm = 1.0 + 0.014 * msq;
    *hs = (*hs + 0.028 * msq) / fm;
    *hs_hk = (*hs_hk) / fm;
    *hs_rt = (*hs_rt) / fm;
    *hs_msq = 0.028 / fm - 0.014 * *hs / fm;
}

void cft(xfoil_t *xf, double hk, double rt, double msq, double *cf, double *cf_hk, double *cf_rt, double *cf_msq)
{
    static const double gam = 1.4;

    /* turbulent skin friction function (Cf), Coles */
    const double gm1 = gam - 1.0;
    const double fc = sqrt(1.0 + 0.5 * gm1 * msq);
    double grt = log(rt / fc);
    grt = fmax_f(grt, 3.0);

    const double gex = -1.74 - 0.31 * hk;

    double arg = -(1.33 * hk);
    arg = fmax_f(-20.0, arg);

    const double thk = tanh(4.0 - hk / 0.875);

    const double cfo = xf->cffac * 0.3 * exp(arg) * pow(grt / 2.3026, gex);
    *cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc;
    *cf_hk = (-(1.33 * cfo) - 0.31 * log(grt / 2.3026) * cfo - 1.1e-4 * (1.0 - P2(thk)) / 0.875) / fc;
    *cf_rt = gex * cfo / (fc * grt) / rt;
    *cf_msq = gex * cfo / (fc * grt) * (-(0.25 * gm1 / P2(fc))) - 0.25 * gm1 * *cf / P2(fc);
}

void hct(double hk, double msq, double *hc, double *hc_hk, double *hc_msq)
{
    /* density shape parameter, from Whitfield */
    *hc = msq * (0.064 / (hk - 0.8) + 0.251);
    *hc_hk = msq * (-(0.064 / P2(hk - 0.8)));
    *hc_msq = 0.064 / (hk - 0.8) + 0.251;
}

/* ------------------------------------------------------------------------- */
/* Part 2: amplification rates (dampl, dampl2, axset) */

void dampl(double hk, double th, double rt, double *ax, double *ax_hk, double *ax_th, double *ax_rt)
{
    static const double dgr = 0.08;

    const double hmi = 1.0 / (hk - 1.0);
    const double hmi_hk = -P2(hmi);

    /* log10(critical Rth) - H correlation for Falkner-Skan profiles */
    const double aa = 2.492 * pow(hmi, 0.43);
    const double aa_hk = (aa / hmi) * 0.43 * hmi_hk;

    const double bb = tanh(14.0 * hmi - 9.24);
    const double bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;

    const double grcrit = aa + 0.7 * (bb + 1.0);
    const double grc_hk = aa_hk + 0.7 * bb_hk;

    const double gr = log10(rt);
    const double gr_rt = 1.0 / (2.3025851 * rt);

    if (gr < grcrit - dgr) {
        /* no amplification for Rtheta < Rcrit */
        *ax = 0.0;
        *ax_hk = 0.0;
        *ax_th = 0.0;
        *ax_rt = 0.0;
    } else {
        /* steep cubic ramp turning AX on smoothly between -DGR < log10(Rtheta/Rcrit) < DGR */
        const double rnorm = (gr - (grcrit - dgr)) / (2.0 * dgr);
        const double rn_hk = -(grc_hk / (2.0 * dgr));
        const double rn_rt = gr_rt / (2.0 * dgr);

        double rfac, rfac_hk, rfac_rt;
        if (rnorm >= 1.0) {
            rfac = 1.0;
            rfac_hk = 0.0;
            rfac_rt = 0.0;
        } else {
            rfac = 3.0 * P2(rnorm) - 2.0 * P3(rnorm);
            const double rfac_rn = 6.0 * rnorm - 6.0 * P2(rnorm);
            rfac_hk = rfac_rn * rn_hk;
            rfac_rt = rfac_rn * rn_rt;
        }

        /* amplification envelope slope correlation for Falkner-Skan */
        const double arg = 3.87 * hmi - 2.52;
        const double arg_hk = 3.87 * hmi_hk;

        const double ex = exp(-P2(arg));
        const double ex_hk = ex * (-(2.0 * arg * arg_hk));

        const double dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
        const double dadr_hk = 0.028 - 0.0345 * ex_hk;

        /* new m(H) correlation, 1 March 91 */
        const double af = -0.05 + 2.7 * hmi - 5.5 * P2(hmi) + 3.0 * P3(hmi);
        const double af_hmi = 2.7 - 11.0 * hmi + 9.0 * P2(hmi);
        const double af_hk = af_hmi * hmi_hk;

        *ax = (af * dadr / th) * rfac;
        *ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
        *ax_th = -(*ax / th);
        *ax_rt = (af * dadr / th) * rfac_rt;
    }
}

void dampl2(double hk, double th, double rt, double *ax, double *ax_hk, double *ax_th, double *ax_rt)
{
    static const double dgr = 0.08;
    static const double hk1 = 3.5, hk2 = 4.0;

    const double hmi = 1.0 / (hk - 1.0);
    const double hmi_hk = -P2(hmi);

    /* log10(critical Rth) - H correlation for Falkner-Skan profiles */
    const double aa = 2.492 * pow(hmi, 0.43);
    const double aa_hk = (aa / hmi) * 0.43 * hmi_hk;

    const double bb = tanh(14.0 * hmi - 9.24);
    const double bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;

    const double grc = aa + 0.7 * (bb + 1.0);
    const double grc_hk = aa_hk + 0.7 * bb_hk;

    const double gr = log10(rt);
    const double gr_rt = 1.0 / (2.3025851 * rt);

    if (gr < grc - dgr) {
        /* no amplification for Rtheta < Rcrit */
        *ax = 0.0;
        *ax_hk = 0.0;
        *ax_th = 0.0;
        *ax_rt = 0.0;
    } else {
        /* steep cubic ramp turning AX on smoothly as Rtheta exceeds Rcrit */
        const double rnorm = (gr - (grc - dgr)) / (2.0 * dgr);
        const double rn_hk = -(grc_hk / (2.0 * dgr));
        const double rn_rt = gr_rt / (2.0 * dgr);

        double rfac, rfac_hk, rfac_rt;
        if (rnorm >= 1.0) {
            rfac = 1.0;
            rfac_hk = 0.0;
            rfac_rt = 0.0;
        } else {
            rfac = 3.0 * P2(rnorm) - 2.0 * P3(rnorm);
            const double rfac_rn = 6.0 * rnorm - 6.0 * P2(rnorm);
            rfac_hk = rfac_rn * rn_hk;
            rfac_rt = rfac_rn * rn_rt;
        }

        /* envelope amplification rate with respect to Rtheta: DADR = dN/dRtheta = f(H) */
        const double arg = 3.87 * hmi - 2.52;
        const double arg_hk = 3.87 * hmi_hk;

        const double ex = exp(-P2(arg));
        const double ex_hk = ex * (-(2.0 * arg * arg_hk));

        const double dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
        const double dadr_hk = 0.028 - 0.0345 * ex_hk;

        /* conversion factor from d/dRtheta to d/dx: AF = Theta dRtheta/dx = f(H) */
        const double brg = -(20.0 * hmi);
        const double af = -0.05 + 2.7 * hmi - 5.5 * P2(hmi) + 3.0 * P3(hmi) + 0.1 * exp(brg);
        const double af_hmi = 2.7 - 11.0 * hmi + 9.0 * P2(hmi) - 2.0 * exp(brg);
        const double af_hk = af_hmi * hmi_hk;

        /* amplification rate with respect to x; RFAC shuts it off below Rcrit */
        *ax = (af * dadr / th) * rfac;
        *ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
        *ax_th = -(*ax / th);
        *ax_rt = (af * dadr / th) * rfac_rt;
    }

    if (hk < hk1) return;

    /* non-envelope max-amplification correction for separated profiles */
    const double hnorm = (hk - hk1) / (hk2 - hk1);
    const double hn_hk = 1.0 / (hk2 - hk1);

    /* blending fraction HFAC = 0..1 over HK1 < HK < HK2 */
    double hfac, hf_hk;
    if (hnorm >= 1.0) {
        hfac = 1.0;
        hf_hk = 0.0;
    } else {
        hfac = 3.0 * P2(hnorm) - 2.0 * P3(hnorm);
        hf_hk = (6.0 * hnorm - 6.0 * P2(hnorm)) * hn_hk;
    }

    /* "normal" envelope amplification rate AX1 */
    const double ax1 = *ax;
    const double ax1_hk = *ax_hk;
    const double ax1_th = *ax_th;
    const double ax1_rt = *ax_rt;

    /* modified amplification rate AX2 */
    const double gr0 = 0.30 + 0.35 * exp(-(0.15 * (hk - 5.0)));
    const double gr0_hk = -(0.35 * exp(-(0.15 * (hk - 5.0))) * 0.15);

    const double tnr = tanh(1.2 * (gr - gr0));
    const double tnr_rt = (1.0 - P2(tnr)) * 1.2 * gr_rt;
    const double tnr_hk = -((1.0 - P2(tnr)) * 1.2 * gr0_hk);

    double ax2 = (0.086 * tnr - 0.25 / pow(hk - 1.0, 1.5)) / th;
    double ax2_hk = (0.086 * tnr_hk + 1.5 * 0.25 / pow(hk - 1.0, 2.5)) / th;
    double ax2_rt = (0.086 * tnr_rt) / th;
    double ax2_th = -(ax2 / th);

    if (ax2 < 0.0) {
        ax2 = 0.0;
        ax2_hk = 0.0;
        ax2_rt = 0.0;
        ax2_th = 0.0;
    }

    /* blend the two amplification rates */
    *ax = hfac * ax2 + (1.0 - hfac) * ax1;
    *ax_hk = hfac * ax2_hk + (1.0 - hfac) * ax1_hk + hf_hk * (ax2 - ax1);
    *ax_rt = hfac * ax2_rt + (1.0 - hfac) * ax1_rt;
    *ax_th = hfac * ax2_th + (1.0 - hfac) * ax1_th;
}

void axset(double hk1, double t1, double rt1, double a1, double hk2, double t2, double rt2, double a2,
           double acrit, int idampv, double *ax, double *ax_hk1, double *ax_t1, double *ax_rt1,
           double *ax_a1, double *ax_hk2, double *ax_t2, double *ax_rt2, double *ax_a2)
{
    double ax1, ax1_hk1, ax1_t1, ax1_rt1, ax2, ax2_hk2, ax2_t2, ax2_rt2;

    /* 2nd-order */
    if (idampv == 0) {
        dampl(hk1, t1, rt1, &ax1, &ax1_hk1, &ax1_t1, &ax1_rt1);
        dampl(hk2, t2, rt2, &ax2, &ax2_hk2, &ax2_t2, &ax2_rt2);
    } else {
        dampl2(hk1, t1, rt1, &ax1, &ax1_hk1, &ax1_t1, &ax1_rt1);
        dampl2(hk2, t2, rt2, &ax2, &ax2_hk2, &ax2_t2, &ax2_rt2);
    }

    /* rms-average version (seems a little better on coarse grids) */
    const double axsq = 0.5 * (P2(ax1) + P2(ax2));
    double axa, axa_ax1, axa_ax2;
    if (axsq <= 0.0) {
        axa = 0.0;
        axa_ax1 = 0.0;
        axa_ax2 = 0.0;
    } else {
        axa = sqrt(axsq);
        axa_ax1 = 0.5 * ax1 / axa;
        axa_ax2 = 0.5 * ax2 / axa;
    }

    /* small additional term to ensure dN/dx > 0 near N = Ncrit */
    const double arg = fmin_f(20.0 * (acrit - 0.5 * (a1 + a2)), 20.0);
    double exn, exn_a1, exn_a2;
    if (arg <= 0.0) {
        exn = 1.0;
        exn_a1 = 0.0;
        exn_a2 = 0.0;
    } else {
        exn = exp(-arg);
        exn_a1 = 20.0 * 0.5 * exn;
        exn_a2 = 20.0 * 0.5 * exn;
    }

    const double dax = exn * 0.002 / (t1 + t2);
    const double dax_a1 = exn_a1 * 0.002 / (t1 + t2);
    const double dax_a2 = exn_a2 * 0.002 / (t1 + t2);
    const double dax_t1 = -(dax / (t1 + t2));
    const double dax_t2 = -(dax / (t1 + t2));

    *ax = axa + dax;

    *ax_hk1 = axa_ax1 * ax1_hk1;
    *ax_t1 = axa_ax1 * ax1_t1 + dax_t1;
    *ax_rt1 = axa_ax1 * ax1_rt1;
    *ax_a1 = dax_a1;

    *ax_hk2 = axa_ax2 * ax2_hk2;
    *ax_t2 = axa_ax2 * ax2_t2 + dax_t2;
    *ax_rt2 = axa_ax2 * ax2_rt2;
    *ax_a2 = dax_a2;
}

/* ------------------------------------------------------------------------- */
/* Part 3: station variables. Fortran "2" quantities are com2.<field> and "1"
 * quantities com1.<field> (X2 -> com2.x, HK2_U2 -> com2.hk_u, DI2_S2 -> com2.di_s). */

void blprv(xfoil_t *xf, double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei)
{
    xf->com2.x = xsi;
    xf->com2.ampl = ami;
    xf->com2.s = cti;
    xf->com2.t = thi;
    xf->com2.d = dsi - dswaki;
    xf->com2.dw = dswaki;

    xf->com2.u = uei * (1.0 - xf->tkbl) / (1.0 - xf->tkbl * P2(uei / xf->qinfbl));
    xf->com2.u_uei = (1.0 + xf->tkbl * (2.0 * xf->com2.u * uei / P2(xf->qinfbl) - 1.0)) / (1.0 - xf->tkbl * P2(uei / xf->qinfbl));
    xf->com2.u_ms = (xf->com2.u * P2(uei / xf->qinfbl) - uei) * xf->tkbl_ms / (1.0 - xf->tkbl * P2(uei / xf->qinfbl));
}

void blkin(xfoil_t *xf)
{
    double hk2_h2, hk2_m2;

    /* edge Mach number ** 2 */
    xf->com2.m = xf->com2.u * xf->com2.u * xf->hstinv / (xf->gm1bl * (1.0 - 0.5 * xf->com2.u * xf->com2.u * xf->hstinv));
    const double tr2 = 1.0 + 0.5 * xf->gm1bl * xf->com2.m;
    xf->com2.m_u = 2.0 * xf->com2.m * tr2 / xf->com2.u;
    xf->com2.m_ms = xf->com2.u * xf->com2.u * tr2 / (xf->gm1bl * (1.0 - 0.5 * xf->com2.u * xf->com2.u * xf->hstinv)) * xf->hstinv_ms;

    /* edge static density (isentropic relation) */
    xf->com2.r = xf->rstbl * pow(tr2, -(1.0 / xf->gm1bl));
    xf->com2.r_u = -(xf->com2.r / tr2 * 0.5 * xf->com2.m_u);
    xf->com2.r_ms = -(xf->com2.r / tr2 * 0.5 * xf->com2.m_ms) + xf->rstbl_ms * pow(tr2, -(1.0 / xf->gm1bl));

    /* shape parameter */
    xf->com2.h = xf->com2.d / xf->com2.t;
    xf->com2.h_d = 1.0 / xf->com2.t;
    xf->com2.h_t = -(xf->com2.h / xf->com2.t);

    /* edge static/stagnation enthalpy */
    const double herat = 1.0 - 0.5 * xf->com2.u * xf->com2.u * xf->hstinv;
    const double he_u2 = -(xf->com2.u * xf->hstinv);
    const double he_ms = -(0.5 * xf->com2.u * xf->com2.u * xf->hstinv_ms);

    /* molecular viscosity */
    xf->com2.v = sqrt(P3(herat)) * (1.0 + xf->hvrat) / (herat + xf->hvrat) / xf->reybl;
    const double v2_he = xf->com2.v * (1.5 / herat - 1.0 / (herat + xf->hvrat));

    xf->com2.v_u = v2_he * he_u2;
    xf->com2.v_ms = -(xf->com2.v / xf->reybl * xf->reybl_ms) + v2_he * he_ms;
    xf->com2.v_re = -(xf->com2.v / xf->reybl * xf->reybl_re);

    /* kinematic shape parameter */
    hkin(xf->com2.h, xf->com2.m, &xf->com2.hk, &hk2_h2, &hk2_m2);

    xf->com2.hk_u = hk2_m2 * xf->com2.m_u;
    xf->com2.hk_t = hk2_h2 * xf->com2.h_t;
    xf->com2.hk_d = hk2_h2 * xf->com2.h_d;
    xf->com2.hk_ms = hk2_m2 * xf->com2.m_ms;

    /* momentum thickness Reynolds number */
    xf->com2.rt = xf->com2.r * xf->com2.u * xf->com2.t / xf->com2.v;
    xf->com2.rt_u = xf->com2.rt * (1.0 / xf->com2.u + xf->com2.r_u / xf->com2.r - xf->com2.v_u / xf->com2.v);
    xf->com2.rt_t = xf->com2.rt / xf->com2.t;
    xf->com2.rt_ms = xf->com2.rt * (xf->com2.r_ms / xf->com2.r - xf->com2.v_ms / xf->com2.v);
    xf->com2.rt_re = xf->com2.rt * (-(xf->com2.v_re / xf->com2.v));
}

void blvar(xfoil_t *xf, int ityp)
{
    double hc2_hk2, hc2_m2, hs2_hk2, hs2_rt2, hs2_m2;
    double cf2_hk2, cf2_rt2, cf2_m2;

    if (ityp == 3) xf->com2.hk = fmax_f(xf->com2.hk, 1.00005);
    if (ityp != 3) xf->com2.hk = fmax_f(xf->com2.hk, 1.05000);

    /* density thickness shape parameter (H**) */
    hct(xf->com2.hk, xf->com2.m, &xf->com2.hc, &hc2_hk2, &hc2_m2);
    xf->com2.hc_u = hc2_hk2 * xf->com2.hk_u + hc2_m2 * xf->com2.m_u;
    xf->com2.hc_t = hc2_hk2 * xf->com2.hk_t;
    xf->com2.hc_d = hc2_hk2 * xf->com2.hk_d;
    xf->com2.hc_ms = hc2_hk2 * xf->com2.hk_ms + hc2_m2 * xf->com2.m_ms;

    /* KE thickness shape parameter from H - H* correlations */
    if (ityp == 1)
        hsl(xf->com2.hk, xf->com2.rt, xf->com2.m, &xf->com2.hs, &hs2_hk2, &hs2_rt2, &hs2_m2);
    else
        hst(xf->com2.hk, xf->com2.rt, xf->com2.m, &xf->com2.hs, &hs2_hk2, &hs2_rt2, &hs2_m2);

    xf->com2.hs_u = hs2_hk2 * xf->com2.hk_u + hs2_rt2 * xf->com2.rt_u + hs2_m2 * xf->com2.m_u;
    xf->com2.hs_t = hs2_hk2 * xf->com2.hk_t + hs2_rt2 * xf->com2.rt_t;
    xf->com2.hs_d = hs2_hk2 * xf->com2.hk_d;
    xf->com2.hs_ms = hs2_hk2 * xf->com2.hk_ms + hs2_rt2 * xf->com2.rt_ms + hs2_m2 * xf->com2.m_ms;
    xf->com2.hs_re = hs2_rt2 * xf->com2.rt_re;

    /* normalized slip velocity Us */
    xf->com2.us = 0.5 * xf->com2.hs * (1.0 - (xf->com2.hk - 1.0) / (xf->gbcon * xf->com2.h));
    const double us2_hs2 = 0.5 * (1.0 - (xf->com2.hk - 1.0) / (xf->gbcon * xf->com2.h));
    const double us2_hk2 = 0.5 * xf->com2.hs * (-(1.0 / (xf->gbcon * xf->com2.h)));
    const double us2_h2 = 0.5 * xf->com2.hs * (xf->com2.hk - 1.0) / (xf->gbcon * P2(xf->com2.h));

    xf->com2.us_u = us2_hs2 * xf->com2.hs_u + us2_hk2 * xf->com2.hk_u;
    xf->com2.us_t = us2_hs2 * xf->com2.hs_t + us2_hk2 * xf->com2.hk_t + us2_h2 * xf->com2.h_t;
    xf->com2.us_d = us2_hs2 * xf->com2.hs_d + us2_hk2 * xf->com2.hk_d + us2_h2 * xf->com2.h_d;
    xf->com2.us_ms = us2_hs2 * xf->com2.hs_ms + us2_hk2 * xf->com2.hk_ms;
    xf->com2.us_re = us2_hs2 * xf->com2.hs_re;

    if (ityp <= 2 && xf->com2.us > 0.95) {
        xf->com2.us = 0.98;
        xf->com2.us_u = 0.0;
        xf->com2.us_t = 0.0;
        xf->com2.us_d = 0.0;
        xf->com2.us_ms = 0.0;
        xf->com2.us_re = 0.0;
    }

    if (ityp == 3 && xf->com2.us > 0.99995) {
        xf->com2.us = 0.99995;
        xf->com2.us_u = 0.0;
        xf->com2.us_t = 0.0;
        xf->com2.us_d = 0.0;
        xf->com2.us_ms = 0.0;
        xf->com2.us_re = 0.0;
    }

    /* equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2  (new 12 Oct 94) */
    double gcc = 0.0;
    double hkc = xf->com2.hk - 1.0;
    double hkc_hk2 = 1.0;
    double hkc_rt2 = 0.0;
    if (ityp == 2) {
        gcc = xf->gccon;
        hkc = xf->com2.hk - 1.0 - gcc / xf->com2.rt;
        hkc_hk2 = 1.0;
        hkc_rt2 = gcc / P2(xf->com2.rt);
        if (hkc < 0.01) {
            hkc = 0.01;
            hkc_hk2 = 0.0;
            hkc_rt2 = 0.0;
        }
    }

    const double hkb = xf->com2.hk - 1.0;
    const double usb = 1.0 - xf->com2.us;
    xf->com2.cq = sqrt(xf->ctcon * xf->com2.hs * hkb * P2(hkc) / (usb * xf->com2.h * P2(xf->com2.hk)));
    const double cq2_hs2 = xf->ctcon * hkb * P2(hkc) / (usb * xf->com2.h * P2(xf->com2.hk)) * 0.5 / xf->com2.cq;
    const double cq2_us2 = xf->ctcon * xf->com2.hs * hkb * P2(hkc) / (usb * xf->com2.h * P2(xf->com2.hk)) / usb * 0.5 / xf->com2.cq;
    const double cq2_hk2 = xf->ctcon * xf->com2.hs * P2(hkc) / (usb * xf->com2.h * P2(xf->com2.hk)) * 0.5 / xf->com2.cq
                           - xf->ctcon * xf->com2.hs * hkb * P2(hkc) / (usb * xf->com2.h * P3(xf->com2.hk)) * 2.0 * 0.5 / xf->com2.cq
                           + xf->ctcon * xf->com2.hs * hkb * hkc / (usb * xf->com2.h * P2(xf->com2.hk)) * 2.0 * 0.5 / xf->com2.cq * hkc_hk2;
    const double cq2_rt2 = xf->ctcon * xf->com2.hs * hkb * hkc / (usb * xf->com2.h * P2(xf->com2.hk)) * 2.0 * 0.5 / xf->com2.cq * hkc_rt2;
    const double cq2_h2 = -(xf->ctcon * xf->com2.hs * hkb * P2(hkc) / (usb * xf->com2.h * P2(xf->com2.hk)) / xf->com2.h * 0.5 / xf->com2.cq);

    xf->com2.cq_u = cq2_hs2 * xf->com2.hs_u + cq2_us2 * xf->com2.us_u + cq2_hk2 * xf->com2.hk_u;
    xf->com2.cq_t = cq2_hs2 * xf->com2.hs_t + cq2_us2 * xf->com2.us_t + cq2_hk2 * xf->com2.hk_t;
    xf->com2.cq_d = cq2_hs2 * xf->com2.hs_d + cq2_us2 * xf->com2.us_d + cq2_hk2 * xf->com2.hk_d;
    xf->com2.cq_ms = cq2_hs2 * xf->com2.hs_ms + cq2_us2 * xf->com2.us_ms + cq2_hk2 * xf->com2.hk_ms;
    xf->com2.cq_re = cq2_hs2 * xf->com2.hs_re + cq2_us2 * xf->com2.us_re;

    xf->com2.cq_u = xf->com2.cq_u + cq2_rt2 * xf->com2.rt_u;
    xf->com2.cq_t = xf->com2.cq_t + cq2_h2 * xf->com2.h_t + cq2_rt2 * xf->com2.rt_t;
    xf->com2.cq_d = xf->com2.cq_d + cq2_h2 * xf->com2.h_d;
    xf->com2.cq_ms = xf->com2.cq_ms + cq2_rt2 * xf->com2.rt_ms;
    xf->com2.cq_re = xf->com2.cq_re + cq2_rt2 * xf->com2.rt_re;

    /* skin friction coefficient */
    if (ityp == 3) {
        /* wake */
        xf->com2.cf = 0.0;
        cf2_hk2 = 0.0;
        cf2_rt2 = 0.0;
        cf2_m2 = 0.0;
    } else if (ityp == 1) {
        /* laminar */
        cfl(xf->com2.hk, xf->com2.rt, xf->com2.m, &xf->com2.cf, &cf2_hk2, &cf2_rt2, &cf2_m2);
    } else {
        /* turbulent; use laminar Cf if it is larger (only for unreasonably small Rtheta) */
        double cf2l, cf2l_hk2, cf2l_rt2, cf2l_m2;
        cft(xf, xf->com2.hk, xf->com2.rt, xf->com2.m, &xf->com2.cf, &cf2_hk2, &cf2_rt2, &cf2_m2);
        cfl(xf->com2.hk, xf->com2.rt, xf->com2.m, &cf2l, &cf2l_hk2, &cf2l_rt2, &cf2l_m2);
        if (cf2l > xf->com2.cf) {
            xf->com2.cf = cf2l;
            cf2_hk2 = cf2l_hk2;
            cf2_rt2 = cf2l_rt2;
            cf2_m2 = cf2l_m2;
        }
    }

    xf->com2.cf_u = cf2_hk2 * xf->com2.hk_u + cf2_rt2 * xf->com2.rt_u + cf2_m2 * xf->com2.m_u;
    xf->com2.cf_t = cf2_hk2 * xf->com2.hk_t + cf2_rt2 * xf->com2.rt_t;
    xf->com2.cf_d = cf2_hk2 * xf->com2.hk_d;
    xf->com2.cf_ms = cf2_hk2 * xf->com2.hk_ms + cf2_rt2 * xf->com2.rt_ms + cf2_m2 * xf->com2.m_ms;
    xf->com2.cf_re = cf2_rt2 * xf->com2.rt_re;

    /* dissipation function 2 CD / H* */
    if (ityp == 1) {
        /* laminar */
        double di2_hk2, di2_rt2;
        dil(xf->com2.hk, xf->com2.rt, &xf->com2.di, &di2_hk2, &di2_rt2);

        xf->com2.di_u = di2_hk2 * xf->com2.hk_u + di2_rt2 * xf->com2.rt_u;
        xf->com2.di_t = di2_hk2 * xf->com2.hk_t + di2_rt2 * xf->com2.rt_t;
        xf->com2.di_d = di2_hk2 * xf->com2.hk_d;
        xf->com2.di_s = 0.0;
        xf->com2.di_ms = di2_hk2 * xf->com2.hk_ms + di2_rt2 * xf->com2.rt_ms;
        xf->com2.di_re = di2_rt2 * xf->com2.rt_re;
    } else if (ityp == 2) {
        /* turbulent wall contribution */
        double cf2t, cf2t_hk2, cf2t_rt2, cf2t_m2;
        cft(xf, xf->com2.hk, xf->com2.rt, xf->com2.m, &cf2t, &cf2t_hk2, &cf2t_rt2, &cf2t_m2);
        const double cf2t_u2 = cf2t_hk2 * xf->com2.hk_u + cf2t_rt2 * xf->com2.rt_u + cf2t_m2 * xf->com2.m_u;
        const double cf2t_t2 = cf2t_hk2 * xf->com2.hk_t + cf2t_rt2 * xf->com2.rt_t;
        const double cf2t_d2 = cf2t_hk2 * xf->com2.hk_d;
        const double cf2t_ms = cf2t_hk2 * xf->com2.hk_ms + cf2t_rt2 * xf->com2.rt_ms + cf2t_m2 * xf->com2.m_ms;
        const double cf2t_re = cf2t_rt2 * xf->com2.rt_re;

        xf->com2.di = (0.5 * cf2t * xf->com2.us) * 2.0 / xf->com2.hs;
        const double di2_hs2 = -((0.5 * cf2t * xf->com2.us) * 2.0 / P2(xf->com2.hs));
        const double di2_us2 = (0.5 * cf2t) * 2.0 / xf->com2.hs;
        const double di2_cf2t = (0.5 * xf->com2.us) * 2.0 / xf->com2.hs;

        xf->com2.di_s = 0.0;
        xf->com2.di_u = di2_hs2 * xf->com2.hs_u + di2_us2 * xf->com2.us_u + di2_cf2t * cf2t_u2;
        xf->com2.di_t = di2_hs2 * xf->com2.hs_t + di2_us2 * xf->com2.us_t + di2_cf2t * cf2t_t2;
        xf->com2.di_d = di2_hs2 * xf->com2.hs_d + di2_us2 * xf->com2.us_d + di2_cf2t * cf2t_d2;
        xf->com2.di_ms = di2_hs2 * xf->com2.hs_ms + di2_us2 * xf->com2.us_ms + di2_cf2t * cf2t_ms;
        xf->com2.di_re = di2_hs2 * xf->com2.hs_re + di2_us2 * xf->com2.us_re + di2_cf2t * cf2t_re;

        /* minimum Hk for the wake layer to still exist */
        const double grt = log(xf->com2.rt);
        const double hmin = 1.0 + 2.1 / grt;
        const double hm_rt2 = -((2.1 / P2(grt)) / xf->com2.rt);

        /* factor DFAC correcting the wall dissipation for very low Hk */
        const double fl = (xf->com2.hk - 1.0) / (hmin - 1.0);
        const double fl_hk2 = 1.0 / (hmin - 1.0);
        const double fl_rt2 = (-(fl / (hmin - 1.0))) * hm_rt2;

        const double tfl = tanh(fl);
        const double dfac = 0.5 + 0.5 * tfl;
        const double df_fl = 0.5 * (1.0 - P2(tfl));

        const double df_hk2 = df_fl * fl_hk2;
        const double df_rt2 = df_fl * fl_rt2;

        xf->com2.di_s = xf->com2.di_s * dfac;
        xf->com2.di_u = xf->com2.di_u * dfac + xf->com2.di * (df_hk2 * xf->com2.hk_u + df_rt2 * xf->com2.rt_u);
        xf->com2.di_t = xf->com2.di_t * dfac + xf->com2.di * (df_hk2 * xf->com2.hk_t + df_rt2 * xf->com2.rt_t);
        xf->com2.di_d = xf->com2.di_d * dfac + xf->com2.di * (df_hk2 * xf->com2.hk_d);
        xf->com2.di_ms = xf->com2.di_ms * dfac + xf->com2.di * (df_hk2 * xf->com2.hk_ms + df_rt2 * xf->com2.rt_ms);
        xf->com2.di_re = xf->com2.di_re * dfac + xf->com2.di * (df_rt2 * xf->com2.rt_re);
        xf->com2.di = xf->com2.di * dfac;
    } else {
        /* zero wall contribution for the wake */
        xf->com2.di = 0.0;
        xf->com2.di_s = 0.0;
        xf->com2.di_u = 0.0;
        xf->com2.di_t = 0.0;
        xf->com2.di_d = 0.0;
        xf->com2.di_ms = 0.0;
        xf->com2.di_re = 0.0;
    }

    /* turbulent outer layer contribution */
    if (ityp != 1) {
        double dd = P2(xf->com2.s) * (0.995 - xf->com2.us) * 2.0 / xf->com2.hs;
        double dd_hs2 = -(P2(xf->com2.s) * (0.995 - xf->com2.us) * 2.0 / P2(xf->com2.hs));
        double dd_us2 = -(P2(xf->com2.s) * 2.0 / xf->com2.hs);
        const double dd_s2 = xf->com2.s * 2.0 * (0.995 - xf->com2.us) * 2.0 / xf->com2.hs;

        xf->com2.di = xf->com2.di + dd;
        xf->com2.di_s = dd_s2;
        xf->com2.di_u = xf->com2.di_u + dd_hs2 * xf->com2.hs_u + dd_us2 * xf->com2.us_u;
        xf->com2.di_t = xf->com2.di_t + dd_hs2 * xf->com2.hs_t + dd_us2 * xf->com2.us_t;
        xf->com2.di_d = xf->com2.di_d + dd_hs2 * xf->com2.hs_d + dd_us2 * xf->com2.us_d;
        xf->com2.di_ms = xf->com2.di_ms + dd_hs2 * xf->com2.hs_ms + dd_us2 * xf->com2.us_ms;
        xf->com2.di_re = xf->com2.di_re + dd_hs2 * xf->com2.hs_re + dd_us2 * xf->com2.us_re;

        /* laminar stress contribution to the outer layer CD */
        dd = 0.15 * P2(0.995 - xf->com2.us) / xf->com2.rt * 2.0 / xf->com2.hs;
        dd_us2 = -(0.15 * (0.995 - xf->com2.us) * 2.0 / xf->com2.rt * 2.0 / xf->com2.hs);
        dd_hs2 = -(dd / xf->com2.hs);
        const double dd_rt2 = -(dd / xf->com2.rt);

        xf->com2.di = xf->com2.di + dd;
        xf->com2.di_u = xf->com2.di_u + dd_hs2 * xf->com2.hs_u + dd_us2 * xf->com2.us_u + dd_rt2 * xf->com2.rt_u;
        xf->com2.di_t = xf->com2.di_t + dd_hs2 * xf->com2.hs_t + dd_us2 * xf->com2.us_t + dd_rt2 * xf->com2.rt_t;
        xf->com2.di_d = xf->com2.di_d + dd_hs2 * xf->com2.hs_d + dd_us2 * xf->com2.us_d;
        xf->com2.di_ms = xf->com2.di_ms + dd_hs2 * xf->com2.hs_ms + dd_us2 * xf->com2.us_ms + dd_rt2 * xf->com2.rt_ms;
        xf->com2.di_re = xf->com2.di_re + dd_hs2 * xf->com2.hs_re + dd_us2 * xf->com2.us_re + dd_rt2 * xf->com2.rt_re;
    }

    if (ityp == 2) {
        /* use laminar CD if it is larger (only for unreasonably small Rtheta) */
        double di2l, di2l_hk2, di2l_rt2;
        dil(xf->com2.hk, xf->com2.rt, &di2l, &di2l_hk2, &di2l_rt2);
        if (di2l > xf->com2.di) {
            xf->com2.di = di2l;
            xf->com2.di_s = 0.0;
            xf->com2.di_u = di2l_hk2 * xf->com2.hk_u + di2l_rt2 * xf->com2.rt_u;
            xf->com2.di_t = di2l_hk2 * xf->com2.hk_t + di2l_rt2 * xf->com2.rt_t;
            xf->com2.di_d = di2l_hk2 * xf->com2.hk_d;
            xf->com2.di_ms = di2l_hk2 * xf->com2.hk_ms + di2l_rt2 * xf->com2.rt_ms;
            xf->com2.di_re = di2l_rt2 * xf->com2.rt_re;
        }
    }

    if (ityp == 3) {
        /* laminar wake CD, if larger */
        double di2l, di2l_hk2, di2l_rt2;
        dilw(xf->com2.hk, xf->com2.rt, &di2l, &di2l_hk2, &di2l_rt2);
        if (di2l > xf->com2.di) {
            xf->com2.di = di2l;
            xf->com2.di_s = 0.0;
            xf->com2.di_u = di2l_hk2 * xf->com2.hk_u + di2l_rt2 * xf->com2.rt_u;
            xf->com2.di_t = di2l_hk2 * xf->com2.hk_t + di2l_rt2 * xf->com2.rt_t;
            xf->com2.di_d = di2l_hk2 * xf->com2.hk_d;
            xf->com2.di_ms = di2l_hk2 * xf->com2.hk_ms + di2l_rt2 * xf->com2.rt_ms;
            xf->com2.di_re = di2l_rt2 * xf->com2.rt_re;
        }
    }

    if (ityp == 3) {
        /* double dissipation for the wake (two wake halves) */
        xf->com2.di = xf->com2.di * 2.0;
        xf->com2.di_s = xf->com2.di_s * 2.0;
        xf->com2.di_u = xf->com2.di_u * 2.0;
        xf->com2.di_t = xf->com2.di_t * 2.0;
        xf->com2.di_d = xf->com2.di_d * 2.0;
        xf->com2.di_ms = xf->com2.di_ms * 2.0;
        xf->com2.di_re = xf->com2.di_re * 2.0;
    }

    /* BL thickness (Delta) from simplified Green's correlation */
    xf->com2.de = (3.15 + 1.72 / (xf->com2.hk - 1.0)) * xf->com2.t + xf->com2.d;
    const double de2_hk2 = (-(1.72 / P2(xf->com2.hk - 1.0))) * xf->com2.t;

    xf->com2.de_u = de2_hk2 * xf->com2.hk_u;
    xf->com2.de_t = de2_hk2 * xf->com2.hk_t + (3.15 + 1.72 / (xf->com2.hk - 1.0));
    xf->com2.de_d = de2_hk2 * xf->com2.hk_d + 1.0;
    xf->com2.de_ms = de2_hk2 * xf->com2.hk_ms;

    const double hdmax = 12.0;
    if (xf->com2.de > hdmax * xf->com2.t) {
        xf->com2.de = hdmax * xf->com2.t;
        xf->com2.de_u = 0.0;
        xf->com2.de_t = hdmax;
        xf->com2.de_d = 0.0;
        xf->com2.de_ms = 0.0;
    }
}

void blmid(xfoil_t *xf, int ityp)
{
    double cfm_hka, cfm_rta, cfm_ma;

    /* similarity station: "1" variables are the "2" ones */
    if (xf->simi) {
        xf->com1.hk = xf->com2.hk;
        xf->com1.hk_t = xf->com2.hk_t;
        xf->com1.hk_d = xf->com2.hk_d;
        xf->com1.hk_u = xf->com2.hk_u;
        xf->com1.hk_ms = xf->com2.hk_ms;
        xf->com1.rt = xf->com2.rt;
        xf->com1.rt_t = xf->com2.rt_t;
        xf->com1.rt_u = xf->com2.rt_u;
        xf->com1.rt_ms = xf->com2.rt_ms;
        xf->com1.rt_re = xf->com2.rt_re;
        xf->com1.m = xf->com2.m;
        xf->com1.m_u = xf->com2.m_u;
        xf->com1.m_ms = xf->com2.m_ms;
    }

    /* midpoint quantities */
    const double hka = 0.5 * (xf->com1.hk + xf->com2.hk);
    const double rta = 0.5 * (xf->com1.rt + xf->com2.rt);
    const double ma = 0.5 * (xf->com1.m + xf->com2.m);

    /* midpoint skin friction coefficient (zero in the wake) */
    if (ityp == 3) {
        xf->cfm = 0.0;
        cfm_hka = 0.0;
        cfm_rta = 0.0;
        cfm_ma = 0.0;
        xf->cfm_ms = 0.0;
    } else if (ityp == 1) {
        cfl(hka, rta, ma, &xf->cfm, &cfm_hka, &cfm_rta, &cfm_ma);
    } else {
        double cfml, cfml_hka, cfml_rta, cfml_ma;
        cft(xf, hka, rta, ma, &xf->cfm, &cfm_hka, &cfm_rta, &cfm_ma);
        cfl(hka, rta, ma, &cfml, &cfml_hka, &cfml_rta, &cfml_ma);
        if (cfml > xf->cfm) {
            xf->cfm = cfml;
            cfm_hka = cfml_hka;
            cfm_rta = cfml_rta;
            cfm_ma = cfml_ma;
        }
    }

    xf->cfm_u1 = 0.5 * (cfm_hka * xf->com1.hk_u + cfm_ma * xf->com1.m_u + cfm_rta * xf->com1.rt_u);
    xf->cfm_t1 = 0.5 * (cfm_hka * xf->com1.hk_t + cfm_rta * xf->com1.rt_t);
    xf->cfm_d1 = 0.5 * (cfm_hka * xf->com1.hk_d);

    xf->cfm_u2 = 0.5 * (cfm_hka * xf->com2.hk_u + cfm_ma * xf->com2.m_u + cfm_rta * xf->com2.rt_u);
    xf->cfm_t2 = 0.5 * (cfm_hka * xf->com2.hk_t + cfm_rta * xf->com2.rt_t);
    xf->cfm_d2 = 0.5 * (cfm_hka * xf->com2.hk_d);

    xf->cfm_ms = 0.5 * (cfm_hka * xf->com1.hk_ms + cfm_ma * xf->com1.m_ms + cfm_rta * xf->com1.rt_ms
                    + cfm_hka * xf->com2.hk_ms + cfm_ma * xf->com2.m_ms + cfm_rta * xf->com2.rt_ms);
    xf->cfm_re = 0.5 * (cfm_rta * xf->com1.rt_re + cfm_rta * xf->com2.rt_re);
}

/* ------------------------------------------------------------------------- */
/* Part 4: difference equations (bldif, tesys). The Newton-system globals are stored
 * column-major like the Fortran: Fortran VS1(k,l) is vs1[l-1][k-1], i.e. vs1[variable]
 * [equation]; equations 0..3 (amplification/shear lag, momentum, shape, Ue), variables
 * 0..4 (ampl./ctau, theta, delta*, Ue, xi). VSREZ/VSM/VSR/VSX(k) are vsrez/vsm/vsr/vsx[k-1]. */

static void clear_system(xfoil_t *xf)
{
    for (int k = 0; k < 4; k++) {
        xf->vsrez[k] = 0.0;
        xf->vsm[k] = 0.0;
        xf->vsr[k] = 0.0;
        xf->vsx[k] = 0.0;
        for (int l = 0; l < 5; l++) {
            xf->vs1[l][k] = 0.0;
            xf->vs2[l][k] = 0.0;
        }
    }
}

void tesys(xfoil_t *xf, double cte, double tte, double dte)
{
    clear_system(xf);

    blvar(xf, 3);

    xf->vs1[0][0] = -1.0;
    xf->vs2[0][0] = 1.0;
    xf->vsrez[0] = cte - xf->com2.s;

    xf->vs1[1][1] = -1.0;
    xf->vs2[1][1] = 1.0;
    xf->vsrez[1] = tte - xf->com2.t;

    xf->vs1[2][2] = -1.0;
    xf->vs2[2][2] = 1.0;
    xf->vsrez[2] = dte - xf->com2.d - xf->com2.dw;
}

void bldif(xfoil_t *xf, int ityp)
{
    double xlog, ulog, tlog, hlog, ddlog;

    if (ityp == 0) {
        /* similarity logarithmic differences (prescribed) */
        xlog = 1.0;
        ulog = xf->bule;
        tlog = 0.5 * (1.0 - xf->bule);
        hlog = 0.0;
        ddlog = 0.0;
    } else {
        /* usual logarithmic differences */
        xlog = log(xf->com2.x / xf->com1.x);
        ulog = log(xf->com2.u / xf->com1.u);
        tlog = log(xf->com2.t / xf->com1.t);
        hlog = log(xf->com2.hs / xf->com1.hs);
        ddlog = 1.0;
    }

    clear_system(xf);

    /* triggering constant for local upwinding */
    const double hupwt = 1.0;

    double hdcon = 5.0 * hupwt / P2(xf->com2.hk);
    double hd_hk1 = 0.0;
    double hd_hk2 = -(hdcon * 2.0 / xf->com2.hk);

    /* less upwinding in the wake */
    if (ityp == 3) {
        hdcon = hupwt / P2(xf->com2.hk);
        hd_hk1 = 0.0;
        hd_hk2 = -(hdcon * 2.0 / xf->com2.hk);
    }

    /* local upwinding is based on the local change in log(Hk-1) (mainly kicks in at transition) */
    const double arg = fabs((xf->com2.hk - 1.0) / (xf->com1.hk - 1.0));
    const double hl = log(arg);
    const double hl_hk1 = -(1.0 / (xf->com1.hk - 1.0));
    const double hl_hk2 = 1.0 / (xf->com2.hk - 1.0);

    /* local upwinding parameter UPW (0.5 trapezoidal, 1.0 backward Euler) and its linearization */
    const double hlsq = fmin_f(P2(hl), 15.0);
    const double ehh = exp(-(hlsq * hdcon));
    const double upw = 1.0 - 0.5 * ehh;
    const double upw_hl = ehh * hl * hdcon;
    const double upw_hd = 0.5 * ehh * hlsq;

    const double upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
    const double upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

    const double upw_u1 = upw_hk1 * xf->com1.hk_u;
    const double upw_t1 = upw_hk1 * xf->com1.hk_t;
    const double upw_d1 = upw_hk1 * xf->com1.hk_d;
    const double upw_u2 = upw_hk2 * xf->com2.hk_u;
    const double upw_t2 = upw_hk2 * xf->com2.hk_t;
    const double upw_d2 = upw_hk2 * xf->com2.hk_d;
    const double upw_ms = upw_hk1 * xf->com1.hk_ms + upw_hk2 * xf->com2.hk_ms;

    if (ityp == 0) {
        /* LE point: zero amplification factor */
        xf->vs2[0][0] = 1.0;
        xf->vsr[0] = 0.0;
        xf->vsrez[0] = -xf->com2.ampl;
    } else if (ityp == 1) {
        /* laminar part: amplification equation, with the average amplification AX over X1..X2 */
        double ax, ax_hk1, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2;
        axset(xf->com1.hk, xf->com1.t, xf->com1.rt, xf->com1.ampl, xf->com2.hk, xf->com2.t, xf->com2.rt, xf->com2.ampl, xf->amcrit,
              xf->idampv, &ax, &ax_hk1, &ax_t1, &ax_rt1, &ax_a1, &ax_hk2, &ax_t2, &ax_rt2, &ax_a2);

        const double rezc = xf->com2.ampl - xf->com1.ampl - ax * (xf->com2.x - xf->com1.x);
        const double z_ax = -(xf->com2.x - xf->com1.x);

        xf->vs1[0][0] = z_ax * ax_a1 - 1.0;
        xf->vs1[1][0] = z_ax * (ax_hk1 * xf->com1.hk_t + ax_t1 + ax_rt1 * xf->com1.rt_t);
        xf->vs1[2][0] = z_ax * (ax_hk1 * xf->com1.hk_d);
        xf->vs1[3][0] = z_ax * (ax_hk1 * xf->com1.hk_u + ax_rt1 * xf->com1.rt_u);
        xf->vs1[4][0] = ax;
        xf->vs2[0][0] = z_ax * ax_a2 + 1.0;
        xf->vs2[1][0] = z_ax * (ax_hk2 * xf->com2.hk_t + ax_t2 + ax_rt2 * xf->com2.rt_t);
        xf->vs2[2][0] = z_ax * (ax_hk2 * xf->com2.hk_d);
        xf->vs2[3][0] = z_ax * (ax_hk2 * xf->com2.hk_u + ax_rt2 * xf->com2.rt_u);
        xf->vs2[4][0] = -ax;
        xf->vsm[0] = z_ax * (ax_hk1 * xf->com1.hk_ms + ax_rt1 * xf->com1.rt_ms + ax_hk2 * xf->com2.hk_ms + ax_rt2 * xf->com2.rt_ms);
        xf->vsr[0] = z_ax * (ax_rt1 * xf->com1.rt_re + ax_rt2 * xf->com2.rt_re);
        xf->vsx[0] = 0.0;
        xf->vsrez[0] = -rezc;
    } else {
        /* turbulent part: shear lag equation */
        const double sa = (1.0 - upw) * xf->com1.s + upw * xf->com2.s;
        const double cqa = (1.0 - upw) * xf->com1.cq + upw * xf->com2.cq;
        const double cfa = (1.0 - upw) * xf->com1.cf + upw * xf->com2.cf;
        const double hka = (1.0 - upw) * xf->com1.hk + upw * xf->com2.hk;

        const double usa = 0.5 * (xf->com1.us + xf->com2.us);
        const double rta = 0.5 * (xf->com1.rt + xf->com2.rt);
        const double dea = 0.5 * (xf->com1.de + xf->com2.de);
        const double da = 0.5 * (xf->com1.d + xf->com2.d);

        /* increased dissipation length in the wake (decrease its reciprocal) */
        const double ald = (ityp == 3) ? xf->dlcon : 1.0;

        /* equilibrium 1/Ue dUe/dx and its linearization (new 12 Oct 94) */
        double hkc, hkc_hka;
        if (ityp == 2) {
            hkc = hka - 1.0 - xf->gccon / rta;
            hkc_hka = 1.0;
            if (hkc < 0.01) {
                hkc = 0.01;
                hkc_hka = 0.0;
            }
        } else {
            hkc = hka - 1.0;
            hkc_hka = 1.0;
        }

        const double hr = hkc / (xf->gacon * ald * hka);
        const double hr_hka = hkc_hka / (xf->gacon * ald * hka) - hr / hka;

        const double uq = (0.5 * cfa - P2(hr)) / (xf->gbcon * da);
        const double uq_hka = -(2.0 * hr * hr_hka / (xf->gbcon * da));
        const double uq_cfa = 0.5 / (xf->gbcon * da);
        const double uq_da = -(uq / da);
        /* (Fortran also computes the unused sensitivities UQ_T1 .. UQ_RE here; dropped in C6b) */

        const double scc = xf->sccon * 1.333 / (1.0 + usa);
        const double scc_usa = -(scc / (1.0 + usa));

        const double slog = log(xf->com2.s / xf->com1.s);
        const double dxi = xf->com2.x - xf->com1.x;

        const double rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog + dea * 2.0 * (uq * dxi - ulog) * xf->duxcon;

        const double z_cfa = dea * 2.0 * uq_cfa * dxi * xf->duxcon;
        const double z_hka = dea * 2.0 * uq_hka * dxi * xf->duxcon;
        const double z_da = dea * 2.0 * uq_da * dxi * xf->duxcon;
        const double z_sl = -(dea * 2.0);
        const double z_ul = -(dea * 2.0 * xf->duxcon);
        const double z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * xf->duxcon;
        const double z_usa = scc_usa * (cqa - sa * ald) * dxi;
        const double z_cqa = scc * dxi;
        const double z_sa = -(scc * dxi * ald);
        const double z_dea = 2.0 * ((uq * dxi - ulog) * xf->duxcon - slog);

        const double z_upw = z_cqa * (xf->com2.cq - xf->com1.cq) + z_sa * (xf->com2.s - xf->com1.s) + z_cfa * (xf->com2.cf - xf->com1.cf)
                             + z_hka * (xf->com2.hk - xf->com1.hk);
        const double z_de1 = 0.5 * z_dea;
        const double z_de2 = 0.5 * z_dea;
        const double z_us1 = 0.5 * z_usa;
        const double z_us2 = 0.5 * z_usa;
        const double z_d1 = 0.5 * z_da;
        const double z_d2 = 0.5 * z_da;
        const double z_u1 = -(z_ul / xf->com1.u);
        const double z_u2 = z_ul / xf->com2.u;
        const double z_x1 = -z_dxi;
        const double z_x2 = z_dxi;
        const double z_s1 = (1.0 - upw) * z_sa - z_sl / xf->com1.s;
        const double z_s2 = upw * z_sa + z_sl / xf->com2.s;
        const double z_cq1 = (1.0 - upw) * z_cqa;
        const double z_cq2 = upw * z_cqa;
        const double z_cf1 = (1.0 - upw) * z_cfa;
        const double z_cf2 = upw * z_cfa;
        const double z_hk1 = (1.0 - upw) * z_hka;
        const double z_hk2 = upw * z_hka;

        xf->vs1[0][0] = z_s1;
        xf->vs1[1][0] = z_upw * upw_t1 + z_de1 * xf->com1.de_t + z_us1 * xf->com1.us_t;
        xf->vs1[2][0] = z_d1 + z_upw * upw_d1 + z_de1 * xf->com1.de_d + z_us1 * xf->com1.us_d;
        xf->vs1[3][0] = z_u1 + z_upw * upw_u1 + z_de1 * xf->com1.de_u + z_us1 * xf->com1.us_u;
        xf->vs1[4][0] = z_x1;
        xf->vs2[0][0] = z_s2;
        xf->vs2[1][0] = z_upw * upw_t2 + z_de2 * xf->com2.de_t + z_us2 * xf->com2.us_t;
        xf->vs2[2][0] = z_d2 + z_upw * upw_d2 + z_de2 * xf->com2.de_d + z_us2 * xf->com2.us_d;
        xf->vs2[3][0] = z_u2 + z_upw * upw_u2 + z_de2 * xf->com2.de_u + z_us2 * xf->com2.us_u;
        xf->vs2[4][0] = z_x2;
        xf->vsm[0] = z_upw * upw_ms + z_de1 * xf->com1.de_ms + z_us1 * xf->com1.us_ms + z_de2 * xf->com2.de_ms + z_us2 * xf->com2.us_ms;

        xf->vs1[1][0] = xf->vs1[1][0] + z_cq1 * xf->com1.cq_t + z_cf1 * xf->com1.cf_t + z_hk1 * xf->com1.hk_t;
        xf->vs1[2][0] = xf->vs1[2][0] + z_cq1 * xf->com1.cq_d + z_cf1 * xf->com1.cf_d + z_hk1 * xf->com1.hk_d;
        xf->vs1[3][0] = xf->vs1[3][0] + z_cq1 * xf->com1.cq_u + z_cf1 * xf->com1.cf_u + z_hk1 * xf->com1.hk_u;

        xf->vs2[1][0] = xf->vs2[1][0] + z_cq2 * xf->com2.cq_t + z_cf2 * xf->com2.cf_t + z_hk2 * xf->com2.hk_t;
        xf->vs2[2][0] = xf->vs2[2][0] + z_cq2 * xf->com2.cq_d + z_cf2 * xf->com2.cf_d + z_hk2 * xf->com2.hk_d;
        xf->vs2[3][0] = xf->vs2[3][0] + z_cq2 * xf->com2.cq_u + z_cf2 * xf->com2.cf_u + z_hk2 * xf->com2.hk_u;

        xf->vsm[0] = xf->vsm[0] + z_cq1 * xf->com1.cq_ms + z_cf1 * xf->com1.cf_ms
                 + z_hk1 * xf->com1.hk_ms + z_cq2 * xf->com2.cq_ms + z_cf2 * xf->com2.cf_ms + z_hk2 * xf->com2.hk_ms;
        xf->vsr[0] = z_cq1 * xf->com1.cq_re + z_cf1 * xf->com1.cf_re + z_cq2 * xf->com2.cq_re + z_cf2 * xf->com2.cf_re;
        xf->vsx[0] = 0.0;
        xf->vsrez[0] = -rezc;
    }

    /* ---- momentum equation */
    double ha = 0.5 * (xf->com1.h + xf->com2.h);
    const double ma = 0.5 * (xf->com1.m + xf->com2.m);
    const double xa = 0.5 * (xf->com1.x + xf->com2.x);
    const double ta = 0.5 * (xf->com1.t + xf->com2.t);
    double hwa = 0.5 * (xf->com1.dw / xf->com1.t + xf->com2.dw / xf->com2.t);

    /* Cf term, using the central value CFM for better accuracy in drag */
    double cfx = 0.50 * xf->cfm * xa / ta + 0.25 * (xf->com1.cf * xf->com1.x / xf->com1.t + xf->com2.cf * xf->com2.x / xf->com2.t);
    const double cfx_xa = 0.50 * xf->cfm / ta;
    const double cfx_ta = -(0.50 * xf->cfm * xa / P2(ta));

    const double cfx_x1 = 0.25 * xf->com1.cf / xf->com1.t + cfx_xa * 0.5;
    const double cfx_x2 = 0.25 * xf->com2.cf / xf->com2.t + cfx_xa * 0.5;
    const double cfx_t1 = -(0.25 * xf->com1.cf * xf->com1.x / P2(xf->com1.t)) + cfx_ta * 0.5;
    const double cfx_t2 = -(0.25 * xf->com2.cf * xf->com2.x / P2(xf->com2.t)) + cfx_ta * 0.5;
    const double cfx_cf1 = 0.25 * xf->com1.x / xf->com1.t;
    const double cfx_cf2 = 0.25 * xf->com2.x / xf->com2.t;
    const double cfx_cfm = 0.50 * xa / ta;

    double btmp = ha + 2.0 - ma + hwa;

    const double rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;
    double z_cfx = -(xlog * 0.5);
    double z_ha = ulog;
    double z_hwa = ulog;
    const double z_ma = -ulog;
    double z_xl = -(ddlog * 0.5 * cfx);
    double z_ul = ddlog * btmp;
    const double z_tl = ddlog;

    const double z_cfm = z_cfx * cfx_cfm;
    double z_cf1 = z_cfx * cfx_cf1;
    double z_cf2 = z_cfx * cfx_cf2;

    double z_t1 = -(z_tl / xf->com1.t) + z_cfx * cfx_t1 + z_hwa * 0.5 * (-(xf->com1.dw / P2(xf->com1.t)));
    double z_t2 = z_tl / xf->com2.t + z_cfx * cfx_t2 + z_hwa * 0.5 * (-(xf->com2.dw / P2(xf->com2.t)));
    double z_x1 = -(z_xl / xf->com1.x) + z_cfx * cfx_x1;
    double z_x2 = z_xl / xf->com2.x + z_cfx * cfx_x2;
    double z_u1 = -(z_ul / xf->com1.u);
    double z_u2 = z_ul / xf->com2.u;

    xf->vs1[1][1] = 0.5 * z_ha * xf->com1.h_t + z_cfm * xf->cfm_t1 + z_cf1 * xf->com1.cf_t + z_t1;
    xf->vs1[2][1] = 0.5 * z_ha * xf->com1.h_d + z_cfm * xf->cfm_d1 + z_cf1 * xf->com1.cf_d;
    xf->vs1[3][1] = 0.5 * z_ma * xf->com1.m_u + z_cfm * xf->cfm_u1 + z_cf1 * xf->com1.cf_u + z_u1;
    xf->vs1[4][1] = z_x1;
    xf->vs2[1][1] = 0.5 * z_ha * xf->com2.h_t + z_cfm * xf->cfm_t2 + z_cf2 * xf->com2.cf_t + z_t2;
    xf->vs2[2][1] = 0.5 * z_ha * xf->com2.h_d + z_cfm * xf->cfm_d2 + z_cf2 * xf->com2.cf_d;
    xf->vs2[3][1] = 0.5 * z_ma * xf->com2.m_u + z_cfm * xf->cfm_u2 + z_cf2 * xf->com2.cf_u + z_u2;
    xf->vs2[4][1] = z_x2;

    xf->vsm[1] = 0.5 * z_ma * xf->com1.m_ms + z_cfm * xf->cfm_ms + z_cf1 * xf->com1.cf_ms + 0.5 * z_ma * xf->com2.m_ms + z_cf2 * xf->com2.cf_ms;
    xf->vsr[1] = z_cfm * xf->cfm_re + z_cf1 * xf->com1.cf_re + z_cf2 * xf->com2.cf_re;
    xf->vsx[1] = 0.0;
    xf->vsrez[1] = -rezt;

    /* ---- shape parameter equation */
    const double xot1 = xf->com1.x / xf->com1.t;
    const double xot2 = xf->com2.x / xf->com2.t;

    ha = 0.5 * (xf->com1.h + xf->com2.h);
    const double hsa = 0.5 * (xf->com1.hs + xf->com2.hs);
    const double hca = 0.5 * (xf->com1.hc + xf->com2.hc);
    hwa = 0.5 * (xf->com1.dw / xf->com1.t + xf->com2.dw / xf->com2.t);

    const double dix = (1.0 - upw) * xf->com1.di * xot1 + upw * xf->com2.di * xot2;
    cfx = (1.0 - upw) * xf->com1.cf * xot1 + upw * xf->com2.cf * xot2;
    const double dix_upw = xf->com2.di * xot2 - xf->com1.di * xot1;
    const double cfx_upw = xf->com2.cf * xot2 - xf->com1.cf * xot1;

    btmp = 2.0 * hca / hsa + 1.0 - ha - hwa;

    const double rezh = hlog + btmp * ulog + xlog * (0.5 * cfx - dix);
    z_cfx = xlog * 0.5;
    const double z_dix = -xlog;
    const double z_hca = 2.0 * ulog / hsa;
    z_ha = -ulog;
    z_hwa = -ulog;
    z_xl = ddlog * (0.5 * cfx - dix);
    z_ul = ddlog * btmp;
    const double z_hl = ddlog;

    const double z_upw = z_cfx * cfx_upw + z_dix * dix_upw;

    const double z_hs1 = -(hca * ulog / P2(hsa)) - z_hl / xf->com1.hs;
    const double z_hs2 = -(hca * ulog / P2(hsa)) + z_hl / xf->com2.hs;

    z_cf1 = (1.0 - upw) * z_cfx * xot1;
    z_cf2 = upw * z_cfx * xot2;
    const double z_di1 = (1.0 - upw) * z_dix * xot1;
    const double z_di2 = upw * z_dix * xot2;

    z_t1 = (1.0 - upw) * (z_cfx * xf->com1.cf + z_dix * xf->com1.di) * (-(xot1 / xf->com1.t));
    z_t2 = upw * (z_cfx * xf->com2.cf + z_dix * xf->com2.di) * (-(xot2 / xf->com2.t));
    z_x1 = (1.0 - upw) * (z_cfx * xf->com1.cf + z_dix * xf->com1.di) / xf->com1.t - z_xl / xf->com1.x;
    z_x2 = upw * (z_cfx * xf->com2.cf + z_dix * xf->com2.di) / xf->com2.t + z_xl / xf->com2.x;
    z_u1 = -(z_ul / xf->com1.u);
    z_u2 = z_ul / xf->com2.u;

    z_t1 = z_t1 + z_hwa * 0.5 * (-(xf->com1.dw / P2(xf->com1.t)));
    z_t2 = z_t2 + z_hwa * 0.5 * (-(xf->com2.dw / P2(xf->com2.t)));

    xf->vs1[0][2] = z_di1 * xf->com1.di_s;
    xf->vs1[1][2] = z_hs1 * xf->com1.hs_t + z_cf1 * xf->com1.cf_t + z_di1 * xf->com1.di_t + z_t1;
    xf->vs1[2][2] = z_hs1 * xf->com1.hs_d + z_cf1 * xf->com1.cf_d + z_di1 * xf->com1.di_d;
    xf->vs1[3][2] = z_hs1 * xf->com1.hs_u + z_cf1 * xf->com1.cf_u + z_di1 * xf->com1.di_u + z_u1;
    xf->vs1[4][2] = z_x1;
    xf->vs2[0][2] = z_di2 * xf->com2.di_s;
    xf->vs2[1][2] = z_hs2 * xf->com2.hs_t + z_cf2 * xf->com2.cf_t + z_di2 * xf->com2.di_t + z_t2;
    xf->vs2[2][2] = z_hs2 * xf->com2.hs_d + z_cf2 * xf->com2.cf_d + z_di2 * xf->com2.di_d;
    xf->vs2[3][2] = z_hs2 * xf->com2.hs_u + z_cf2 * xf->com2.cf_u + z_di2 * xf->com2.di_u + z_u2;
    xf->vs2[4][2] = z_x2;
    xf->vsm[2] = z_hs1 * xf->com1.hs_ms + z_cf1 * xf->com1.cf_ms + z_di1 * xf->com1.di_ms + z_hs2 * xf->com2.hs_ms + z_cf2 * xf->com2.cf_ms
             + z_di2 * xf->com2.di_ms;
    xf->vsr[2] = z_hs1 * xf->com1.hs_re + z_cf1 * xf->com1.cf_re + z_di1 * xf->com1.di_re + z_hs2 * xf->com2.hs_re + z_cf2 * xf->com2.cf_re
             + z_di2 * xf->com2.di_re;

    xf->vs1[1][2] = xf->vs1[1][2] + 0.5 * (z_hca * xf->com1.hc_t + z_ha * xf->com1.h_t) + z_upw * upw_t1;
    xf->vs1[2][2] = xf->vs1[2][2] + 0.5 * (z_hca * xf->com1.hc_d + z_ha * xf->com1.h_d) + z_upw * upw_d1;
    xf->vs1[3][2] = xf->vs1[3][2] + 0.5 * (z_hca * xf->com1.hc_u) + z_upw * upw_u1;
    xf->vs2[1][2] = xf->vs2[1][2] + 0.5 * (z_hca * xf->com2.hc_t + z_ha * xf->com2.h_t) + z_upw * upw_t2;
    xf->vs2[2][2] = xf->vs2[2][2] + 0.5 * (z_hca * xf->com2.hc_d + z_ha * xf->com2.h_d) + z_upw * upw_d2;
    xf->vs2[3][2] = xf->vs2[3][2] + 0.5 * (z_hca * xf->com2.hc_u) + z_upw * upw_u2;

    xf->vsm[2] = xf->vsm[2] + 0.5 * (z_hca * xf->com1.hc_ms) + z_upw * upw_ms + 0.5 * (z_hca * xf->com2.hc_ms);

    xf->vsx[2] = 0.0;
    xf->vsrez[2] = -rezh;
}

/* ------------------------------------------------------------------------- */
/* Part 5: transition check and system assembly (trchek2, trchek, trdif, blsys) */

int trchek2(xfoil_t *xf)
{
    static const double daeps = 5.0e-5;

    double amplt = 0.0, amplt_a2 = 0.0, ax = 0.0, da2 = 0.0;
    double ax_hk1, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2;
    double ax_hkt = 0.0, ax_tt = 0.0, ax_rtt = 0.0, ax_at = 0.0;
    double wf1 = 0.0, wf1_a1 = 0.0, wf1_a2 = 0.0, wf1_x1 = 0.0, wf1_x2 = 0.0, wf1_xf = 0.0;
    double wf2 = 0.0, wf2_a1 = 0.0, wf2_a2 = 0.0, wf2_x1 = 0.0, wf2_x2 = 0.0, wf2_xf = 0.0;
    double tt = 0.0, dt = 0.0, ut = 0.0, xt_a2 = 0.0, tt_a2 = 0.0, dt_a2 = 0.0, ut_a2 = 0.0;
    double hkt_tt = 0.0, hkt_dt = 0.0, hkt_ut = 0.0, hkt_ms = 0.0;
    double rtt_tt = 0.0, rtt_ut = 0.0, rtt_ms = 0.0, rtt_re = 0.0;

    int result = 1;

    /* save variables and sensitivities at IBL ("2") for future restoration */
    xf->c2sav = xf->com2;

    /* average amplification rate AX over the X1..X2 interval */
    axset(xf->com1.hk, xf->com1.t, xf->com1.rt, xf->com1.ampl, xf->com2.hk, xf->com2.t, xf->com2.rt, xf->com2.ampl, xf->amcrit, xf->idampv,
          &ax, &ax_hk1, &ax_t1, &ax_rt1, &ax_a1, &ax_hk2, &ax_t2, &ax_rt2, &ax_a2);

    /* initial guess for the iterate N2 (AMPL2) at X2 */
    xf->com2.ampl = xf->com1.ampl + ax * (xf->com2.x - xf->com1.x);

    /* solve the implicit system for the amplification AMPL2 */
    int converged = 0;
    for (int itam = 1; itam <= 30; itam++) {
        /* weighting factors WF1, WF2 for defining "T" quantities from 1, 2 */
        double sfa, sfa_a1, sfa_a2;
        if (xf->com2.ampl <= xf->amcrit) {
            /* no transition yet: "T" is the same as "2" */
            amplt = xf->com2.ampl;
            amplt_a2 = 1.0;
            sfa = 1.0;
            sfa_a1 = 0.0;
            sfa_a2 = 0.0;
        } else {
            /* transition in X1..X2: "T" is set from N1, N2 */
            amplt = xf->amcrit;
            amplt_a2 = 0.0;
            sfa = (amplt - xf->com1.ampl) / (xf->com2.ampl - xf->com1.ampl);
            sfa_a1 = (sfa - 1.0) / (xf->com2.ampl - xf->com1.ampl);
            sfa_a2 = (-sfa) / (xf->com2.ampl - xf->com1.ampl);
        }

        double sfx, sfx_x1, sfx_x2, sfx_xf;
        if (xf->xiforc < xf->com2.x) {
            sfx = (xf->xiforc - xf->com1.x) / (xf->com2.x - xf->com1.x);
            sfx_x1 = (sfx - 1.0) / (xf->com2.x - xf->com1.x);
            sfx_x2 = (-sfx) / (xf->com2.x - xf->com1.x);
            sfx_xf = 1.0 / (xf->com2.x - xf->com1.x);
        } else {
            sfx = 1.0;
            sfx_x1 = 0.0;
            sfx_x2 = 0.0;
            sfx_xf = 0.0;
        }

        /* weighting factor from free or forced transition */
        if (sfa < sfx) {
            wf2 = sfa;
            wf2_a1 = sfa_a1;
            wf2_a2 = sfa_a2;
            wf2_x1 = 0.0;
            wf2_x2 = 0.0;
            wf2_xf = 0.0;
        } else {
            wf2 = sfx;
            wf2_a1 = 0.0;
            wf2_a2 = 0.0;
            wf2_x1 = sfx_x1;
            wf2_x2 = sfx_x2;
            wf2_xf = sfx_xf;
        }

        wf1 = 1.0 - wf2;
        wf1_a1 = -wf2_a1;
        wf1_a2 = -wf2_a2;
        wf1_x1 = -wf2_x1;
        wf1_x2 = -wf2_x2;
        wf1_xf = -wf2_xf;

        /* interpolate BL variables to XT */
        xf->xt = xf->com1.x * wf1 + xf->com2.x * wf2;
        tt = xf->com1.t * wf1 + xf->com2.t * wf2;
        dt = xf->com1.d * wf1 + xf->com2.d * wf2;
        ut = xf->com1.u * wf1 + xf->com2.u * wf2;

        xt_a2 = xf->com1.x * wf1_a2 + xf->com2.x * wf2_a2;
        tt_a2 = xf->com1.t * wf1_a2 + xf->com2.t * wf2_a2;
        dt_a2 = xf->com1.d * wf1_a2 + xf->com2.d * wf2_a2;
        ut_a2 = xf->com1.u * wf1_a2 + xf->com2.u * wf2_a2;

        /* temporarily set "2" variables from "T" for BLKIN */
        xf->com2.x = xf->xt;
        xf->com2.t = tt;
        xf->com2.d = dt;
        xf->com2.u = ut;

        /* laminar secondary "T" variables HKT, RTT */
        blkin(xf);

        const double hkt = xf->com2.hk;
        hkt_tt = xf->com2.hk_t;
        hkt_dt = xf->com2.hk_d;
        hkt_ut = xf->com2.hk_u;
        hkt_ms = xf->com2.hk_ms;

        const double rtt = xf->com2.rt;
        rtt_tt = xf->com2.rt_t;
        rtt_ut = xf->com2.rt_u;
        rtt_ms = xf->com2.rt_ms;
        rtt_re = xf->com2.rt_re;

        /* restore clobbered "2" variables, except for AMPL2 */
        const double amsave = xf->com2.ampl;
        xf->com2 = xf->c2sav;
        xf->com2.ampl = amsave;

        /* amplification rate AX over the current X1-XT interval */
        axset(xf->com1.hk, xf->com1.t, xf->com1.rt, xf->com1.ampl, hkt, tt, rtt, amplt, xf->amcrit, xf->idampv,
              &ax, &ax_hk1, &ax_t1, &ax_rt1, &ax_a1, &ax_hkt, &ax_tt, &ax_rtt, &ax_at);

        /* punch out early if there is no amplification here */
        if (ax <= 0.0) {
            converged = 1;
            break;
        }

        /* sensitivity of AX(A2) */
        ax_a2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_a2
                + (ax_hkt * hkt_dt) * dt_a2 + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a2 + ax_at * amplt_a2;

        /* residual for the implicit AMPL2 definition (amplification equation) */
        const double res = xf->com2.ampl - xf->com1.ampl - ax * (xf->com2.x - xf->com1.x);
        const double res_a2 = 1.0 - ax_a2 * (xf->com2.x - xf->com1.x);

        da2 = -(res / res_a2);

        double rlx = 1.0;
        const double dxt = xt_a2 * da2;

        if (rlx * fabs(dxt / (xf->com2.x - xf->com1.x)) > 0.05) rlx = 0.05 * fabs((xf->com2.x - xf->com1.x) / dxt);
        if (rlx * fabs(da2) > 1.0) rlx = 1.0 * fabs(1.0 / da2);

        /* converged? */
        if (fabs(da2) < daeps) {
            converged = 1;
            break;
        }

        if ((xf->com2.ampl > xf->amcrit && xf->com2.ampl + rlx * da2 < xf->amcrit)
            || (xf->com2.ampl < xf->amcrit && xf->com2.ampl + rlx * da2 > xf->amcrit)) {
            /* limited Newton step so AMPL2 doesn't step across AMCRIT either way */
            xf->com2.ampl = xf->amcrit;
        } else {
            /* regular Newton step */
            xf->com2.ampl = xf->com2.ampl + rlx * da2;
        }
    }

    if (!converged) {
        if (xfoil_print) {
            printf(" TRCHEK2: N2 convergence failed.\n");
            printf(" x:%9.5f%9.5f%9.5f  N:%7.3f%7.3f%7.3f  Nx:%8.3f   dN:%10.3E\n", xf->com1.x, xf->xt, xf->com2.x, xf->com1.ampl, amplt, xf->com2.ampl, ax, da2);
        }

        /* if any of the printed variables is NaN, convergence will never be reached */
        if (isnan(xf->com1.x) || isnan(xf->xt) || isnan(xf->com2.x) || isnan(xf->com1.ampl) || isnan(amplt) || isnan(xf->com2.ampl)
            || isnan(ax) || isnan(da2)) {
            result = 0;
            if (ABORT_ON_NAN) return result;
        }
    }
    /* free or forced transition? */
    xf->trfree = xf->com2.ampl >= xf->amcrit;
    xf->trforc = xf->xiforc > xf->com1.x && xf->xiforc <= xf->com2.x;

    /* transition interval flag */
    xf->tran = xf->trforc || xf->trfree;

    if (ABORT_ON_NAN && !xf->tran) return result;

    /* resolve if both forced and free transition */
    if (xf->trfree && xf->trforc) {
        xf->trforc = xf->xiforc < xf->xt;
        xf->trfree = xf->xiforc >= xf->xt;
    }

    if (xf->trforc) {
        /* forced transition: XT is prescribed, the sensitivities are known */
        xf->xt = xf->xiforc;
        xf->xt_a1 = 0.0;
        xf->xt_x1 = 0.0;
        xf->xt_t1 = 0.0;
        xf->xt_d1 = 0.0;
        xf->xt_u1 = 0.0;
        xf->xt_x2 = 0.0;
        xf->xt_t2 = 0.0;
        xf->xt_d2 = 0.0;
        xf->xt_u2 = 0.0;
        xf->xt_ms = 0.0;
        xf->xt_re = 0.0;
        xf->xt_xf = 1.0;
        return result;
    }

    /* free transition: sensitivities of XT, TT, DT, UT */
    xf->xt_x1 = wf1;
    double tt_t1 = wf1;
    double dt_d1 = wf1;
    double ut_u1 = wf1;

    xf->xt_x2 = wf2;
    const double tt_t2 = wf2;
    const double dt_d2 = wf2;
    const double ut_u2 = wf2;

    xf->xt_a1 = xf->com1.x * wf1_a1 + xf->com2.x * wf2_a1;
    const double tt_a1 = xf->com1.t * wf1_a1 + xf->com2.t * wf2_a1;
    const double dt_a1 = xf->com1.d * wf1_a1 + xf->com2.d * wf2_a1;
    const double ut_a1 = xf->com1.u * wf1_a1 + xf->com2.u * wf2_a1;

    xf->xt_x1 = xf->com1.x * wf1_x1 + xf->com2.x * wf2_x1 + xf->xt_x1;
    const double tt_x1 = xf->com1.t * wf1_x1 + xf->com2.t * wf2_x1;
    const double dt_x1 = xf->com1.d * wf1_x1 + xf->com2.d * wf2_x1;
    const double ut_x1 = xf->com1.u * wf1_x1 + xf->com2.u * wf2_x1;

    xf->xt_x2 = xf->com1.x * wf1_x2 + xf->com2.x * wf2_x2 + xf->xt_x2;
    const double tt_x2 = xf->com1.t * wf1_x2 + xf->com2.t * wf2_x2;
    const double dt_x2 = xf->com1.d * wf1_x2 + xf->com2.d * wf2_x2;
    const double ut_x2 = xf->com1.u * wf1_x2 + xf->com2.u * wf2_x2;

    xf->xt_xf = xf->com1.x * wf1_xf + xf->com2.x * wf2_xf;

    /* AX = AX(HK1, T1, RT1, A1, HKT, TT, RTT, AT): sensitivities wrt T1 D1 U1 A1 T2 D2 U2 A2 MS RE */
    ax_t1 = ax_hk1 * xf->com1.hk_t + ax_t1 + ax_rt1 * xf->com1.rt_t + (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_t1;
    const double ax_d1 = ax_hk1 * xf->com1.hk_d + (ax_hkt * hkt_dt) * dt_d1;
    const double ax_u1 = ax_hk1 * xf->com1.hk_u + ax_rt1 * xf->com1.rt_u + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_u1;
    ax_a1 = ax_a1 + (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_a1 + (ax_hkt * hkt_dt) * dt_a1
            + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a1;
    const double ax_x1 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_x1 + (ax_hkt * hkt_dt) * dt_x1
                         + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_x1;

    ax_t2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_t2;
    const double ax_d2 = (ax_hkt * hkt_dt) * dt_d2;
    const double ax_u2 = (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_u2;
    ax_a2 = ax_at * amplt_a2 + (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_a2 + (ax_hkt * hkt_dt) * dt_a2
            + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a2;
    const double ax_x2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_x2 + (ax_hkt * hkt_dt) * dt_x2
                         + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_x2;


    const double ax_ms = ax_hkt * hkt_ms + ax_rtt * rtt_ms + ax_hk1 * xf->com1.hk_ms + ax_rt1 * xf->com1.rt_ms;
    const double ax_re = ax_rtt * rtt_re + ax_rt1 * xf->com1.rt_re;

    /* sensitivities of the residual RES = AMPL2 - AMPL1 - AX*(X2-X1) */
    const double z_ax = -(xf->com2.x - xf->com1.x);

    const double z_a1 = z_ax * ax_a1 - 1.0;
    const double z_t1 = z_ax * ax_t1;
    const double z_d1 = z_ax * ax_d1;
    const double z_u1 = z_ax * ax_u1;
    const double z_x1 = z_ax * ax_x1 + ax;

    const double z_a2 = z_ax * ax_a2 + 1.0;
    const double z_t2 = z_ax * ax_t2;
    const double z_d2 = z_ax * ax_d2;
    const double z_u2 = z_ax * ax_u2;
    const double z_x2 = z_ax * ax_x2 - ax;

    const double z_ms = z_ax * ax_ms;
    const double z_re = z_ax * ax_re;

    /* sensitivities of XT, with RES stationary for the A2 constraint */
    xf->xt_a1 = xf->xt_a1 - (xt_a2 / z_a2) * z_a1;
    xf->xt_t1 = -((xt_a2 / z_a2) * z_t1);
    xf->xt_d1 = -((xt_a2 / z_a2) * z_d1);
    xf->xt_u1 = -((xt_a2 / z_a2) * z_u1);
    xf->xt_x1 = xf->xt_x1 - (xt_a2 / z_a2) * z_x1;
    xf->xt_t2 = -((xt_a2 / z_a2) * z_t2);
    xf->xt_d2 = -((xt_a2 / z_a2) * z_d2);
    xf->xt_u2 = -((xt_a2 / z_a2) * z_u2);
    xf->xt_x2 = xf->xt_x2 - (xt_a2 / z_a2) * z_x2;
    xf->xt_ms = -((xt_a2 / z_a2) * z_ms);
    xf->xt_re = -((xt_a2 / z_a2) * z_re);
    xf->xt_xf = 0.0;

    return result;
}

int trchek(xfoil_t *xf)
{
    /* 2nd-order amplification equation */
    return trchek2(xf);
}

void trdif(xfoil_t *xf)
{
    /* laminar (bl*) and turbulent (bt*) parts of equations 0..2, [equation][variable] */
    double bl1[3][5], bl2[3][5], bt1[3][5], bt2[3][5];
    double blm[3], blr[3], blrez[3], blx[3], btm[3], btr[3], btrez[3], btx[3];

    /* save variables and sensitivities for future restoration */
    xf->c1sav = xf->com1;
    xf->c2sav = xf->com2;

    /* weighting factors for linear interpolation to the transition point */
    const double wf2 = (xf->xt - xf->com1.x) / (xf->com2.x - xf->com1.x);
    const double wf2_xt = 1.0 / (xf->com2.x - xf->com1.x);

    const double wf2_a1 = wf2_xt * xf->xt_a1;
    const double wf2_x1 = wf2_xt * xf->xt_x1 + (wf2 - 1.0) / (xf->com2.x - xf->com1.x);
    const double wf2_x2 = wf2_xt * xf->xt_x2 - wf2 / (xf->com2.x - xf->com1.x);
    const double wf2_t1 = wf2_xt * xf->xt_t1;
    const double wf2_t2 = wf2_xt * xf->xt_t2;
    const double wf2_d1 = wf2_xt * xf->xt_d1;
    const double wf2_d2 = wf2_xt * xf->xt_d2;
    const double wf2_u1 = wf2_xt * xf->xt_u1;
    const double wf2_u2 = wf2_xt * xf->xt_u2;
    const double wf2_ms = wf2_xt * xf->xt_ms;
    const double wf2_re = wf2_xt * xf->xt_re;
    const double wf2_xf = wf2_xt * xf->xt_xf;

    const double wf1 = 1.0 - wf2;
    const double wf1_a1 = -wf2_a1;
    const double wf1_x1 = -wf2_x1;
    const double wf1_x2 = -wf2_x2;
    const double wf1_t1 = -wf2_t1;
    const double wf1_t2 = -wf2_t2;
    const double wf1_d1 = -wf2_d1;
    const double wf1_d2 = -wf2_d2;
    const double wf1_u1 = -wf2_u1;
    const double wf1_u2 = -wf2_u2;
    const double wf1_ms = -wf2_ms;
    const double wf1_re = -wf2_re;
    const double wf1_xf = -wf2_xf;

    /**** FIRST, the laminar part between X1 and XT */

    /* interpolate primary variables to the transition point */
    const double t1 = xf->com1.t, t2 = xf->com2.t, d1 = xf->com1.d, d2 = xf->com2.d, u1 = xf->com1.u, u2 = xf->com2.u;
    const double tt = t1 * wf1 + t2 * wf2;
    const double tt_a1 = t1 * wf1_a1 + t2 * wf2_a1;
    const double tt_x1 = t1 * wf1_x1 + t2 * wf2_x1;
    const double tt_x2 = t1 * wf1_x2 + t2 * wf2_x2;
    const double tt_t1 = t1 * wf1_t1 + t2 * wf2_t1 + wf1;
    const double tt_t2 = t1 * wf1_t2 + t2 * wf2_t2 + wf2;
    const double tt_d1 = t1 * wf1_d1 + t2 * wf2_d1;
    const double tt_d2 = t1 * wf1_d2 + t2 * wf2_d2;
    const double tt_u1 = t1 * wf1_u1 + t2 * wf2_u1;
    const double tt_u2 = t1 * wf1_u2 + t2 * wf2_u2;
    const double tt_ms = t1 * wf1_ms + t2 * wf2_ms;
    const double tt_re = t1 * wf1_re + t2 * wf2_re;
    const double tt_xf = t1 * wf1_xf + t2 * wf2_xf;

    const double dt = d1 * wf1 + d2 * wf2;
    const double dt_a1 = d1 * wf1_a1 + d2 * wf2_a1;
    const double dt_x1 = d1 * wf1_x1 + d2 * wf2_x1;
    const double dt_x2 = d1 * wf1_x2 + d2 * wf2_x2;
    const double dt_t1 = d1 * wf1_t1 + d2 * wf2_t1;
    const double dt_t2 = d1 * wf1_t2 + d2 * wf2_t2;
    const double dt_d1 = d1 * wf1_d1 + d2 * wf2_d1 + wf1;
    const double dt_d2 = d1 * wf1_d2 + d2 * wf2_d2 + wf2;
    const double dt_u1 = d1 * wf1_u1 + d2 * wf2_u1;
    const double dt_u2 = d1 * wf1_u2 + d2 * wf2_u2;
    const double dt_ms = d1 * wf1_ms + d2 * wf2_ms;
    const double dt_re = d1 * wf1_re + d2 * wf2_re;
    const double dt_xf = d1 * wf1_xf + d2 * wf2_xf;

    const double ut = u1 * wf1 + u2 * wf2;
    const double ut_a1 = u1 * wf1_a1 + u2 * wf2_a1;
    const double ut_x1 = u1 * wf1_x1 + u2 * wf2_x1;
    const double ut_x2 = u1 * wf1_x2 + u2 * wf2_x2;
    const double ut_t1 = u1 * wf1_t1 + u2 * wf2_t1;
    const double ut_t2 = u1 * wf1_t2 + u2 * wf2_t2;
    const double ut_d1 = u1 * wf1_d1 + u2 * wf2_d1;
    const double ut_d2 = u1 * wf1_d2 + u2 * wf2_d2;
    const double ut_u1 = u1 * wf1_u1 + u2 * wf2_u1 + wf1;
    const double ut_u2 = u1 * wf1_u2 + u2 * wf2_u2 + wf2;
    const double ut_ms = u1 * wf1_ms + u2 * wf2_ms;
    const double ut_re = u1 * wf1_re + u2 * wf2_re;
    const double ut_xf = u1 * wf1_xf + u2 * wf2_xf;

    /* primary "T" variables at XT (placed into the "2" variables) */
    xf->com2.x = xf->xt;
    xf->com2.t = tt;
    xf->com2.d = dt;
    xf->com2.u = ut;

    xf->com2.ampl = xf->amcrit;
    xf->com2.s = 0.0;

    /* laminar secondary "T" variables, X1-XT midpoint CFM, and the Newton system at X1, XT */
    blkin(xf);
    blvar(xf, 1);
    blmid(xf, 1);
    bldif(xf, 1);

    /* convert sensitivities wrt "T" into sensitivities wrt "1" and "2" variables
     * (the K=1 amplification row is not needed here) */
    for (int k = 1; k < 3; k++) {
        blrez[k] = xf->vsrez[k];
        blm[k] = xf->vsm[k] + xf->vs2[1][k] * tt_ms + xf->vs2[2][k] * dt_ms + xf->vs2[3][k] * ut_ms + xf->vs2[4][k] * xf->xt_ms;
        blr[k] = xf->vsr[k] + xf->vs2[1][k] * tt_re + xf->vs2[2][k] * dt_re + xf->vs2[3][k] * ut_re + xf->vs2[4][k] * xf->xt_re;
        blx[k] = xf->vsx[k] + xf->vs2[1][k] * tt_xf + xf->vs2[2][k] * dt_xf + xf->vs2[3][k] * ut_xf + xf->vs2[4][k] * xf->xt_xf;

        bl1[k][0] = xf->vs1[0][k] + xf->vs2[1][k] * tt_a1 + xf->vs2[2][k] * dt_a1 + xf->vs2[3][k] * ut_a1 + xf->vs2[4][k] * xf->xt_a1;
        bl1[k][1] = xf->vs1[1][k] + xf->vs2[1][k] * tt_t1 + xf->vs2[2][k] * dt_t1 + xf->vs2[3][k] * ut_t1 + xf->vs2[4][k] * xf->xt_t1;
        bl1[k][2] = xf->vs1[2][k] + xf->vs2[1][k] * tt_d1 + xf->vs2[2][k] * dt_d1 + xf->vs2[3][k] * ut_d1 + xf->vs2[4][k] * xf->xt_d1;
        bl1[k][3] = xf->vs1[3][k] + xf->vs2[1][k] * tt_u1 + xf->vs2[2][k] * dt_u1 + xf->vs2[3][k] * ut_u1 + xf->vs2[4][k] * xf->xt_u1;
        bl1[k][4] = xf->vs1[4][k] + xf->vs2[1][k] * tt_x1 + xf->vs2[2][k] * dt_x1 + xf->vs2[3][k] * ut_x1 + xf->vs2[4][k] * xf->xt_x1;

        bl2[k][0] = 0.0;
        bl2[k][1] = xf->vs2[1][k] * tt_t2 + xf->vs2[2][k] * dt_t2 + xf->vs2[3][k] * ut_t2 + xf->vs2[4][k] * xf->xt_t2;
        bl2[k][2] = xf->vs2[1][k] * tt_d2 + xf->vs2[2][k] * dt_d2 + xf->vs2[3][k] * ut_d2 + xf->vs2[4][k] * xf->xt_d2;
        bl2[k][3] = xf->vs2[1][k] * tt_u2 + xf->vs2[2][k] * dt_u2 + xf->vs2[3][k] * ut_u2 + xf->vs2[4][k] * xf->xt_u2;
        bl2[k][4] = xf->vs2[1][k] * tt_x2 + xf->vs2[2][k] * dt_x2 + xf->vs2[3][k] * ut_x2 + xf->vs2[4][k] * xf->xt_x2;
    }

    /**** SECOND, the turbulent part between XT and X2 */

    /* equilibrium shear coefficient CQT at the transition point */
    blvar(xf, 2);

    /* initial shear coefficient value ST at the transition point
     * (CQ2, CQ2_T2, ... are really "CQT", "CQT_TT", ...) */
    const double ctr = xf->ctrcon * exp(-(xf->ctrcex / (xf->com2.hk - 1.0)));
    const double ctr_hk2 = ctr * xf->ctrcex / P2(xf->com2.hk - 1.0);

    const double st = ctr * xf->com2.cq;
    const double st_tt = ctr * xf->com2.cq_t + xf->com2.cq * ctr_hk2 * xf->com2.hk_t;
    const double st_dt = ctr * xf->com2.cq_d + xf->com2.cq * ctr_hk2 * xf->com2.hk_d;
    const double st_ut = ctr * xf->com2.cq_u + xf->com2.cq * ctr_hk2 * xf->com2.hk_u;
    double st_ms = ctr * xf->com2.cq_ms + xf->com2.cq * ctr_hk2 * xf->com2.hk_ms;
    double st_re = ctr * xf->com2.cq_re;

    /* ST sensitivities wrt the actual "1" and "2" variables */
    const double st_a1 = st_tt * tt_a1 + st_dt * dt_a1 + st_ut * ut_a1;
    const double st_x1 = st_tt * tt_x1 + st_dt * dt_x1 + st_ut * ut_x1;
    const double st_x2 = st_tt * tt_x2 + st_dt * dt_x2 + st_ut * ut_x2;
    const double st_t1 = st_tt * tt_t1 + st_dt * dt_t1 + st_ut * ut_t1;
    const double st_t2 = st_tt * tt_t2 + st_dt * dt_t2 + st_ut * ut_t2;
    const double st_d1 = st_tt * tt_d1 + st_dt * dt_d1 + st_ut * ut_d1;
    const double st_d2 = st_tt * tt_d2 + st_dt * dt_d2 + st_ut * ut_d2;
    const double st_u1 = st_tt * tt_u1 + st_dt * dt_u1 + st_ut * ut_u1;
    const double st_u2 = st_tt * tt_u2 + st_dt * dt_u2 + st_ut * ut_u2;
    st_ms = st_tt * tt_ms + st_dt * dt_ms + st_ut * ut_ms + st_ms;
    st_re = st_tt * tt_re + st_dt * dt_re + st_ut * ut_re + st_re;
    const double st_xf = st_tt * tt_xf + st_dt * dt_xf + st_ut * ut_xf;

    xf->com2.ampl = 0.0;
    xf->com2.s = st;

    /* recalculate the turbulent secondary "T" variables using the proper CTI */
    blvar(xf, 2);

    /* "1" variables := "T" variables; reset "2" variables to their saved turbulent values */
    xf->com1 = xf->com2;
    xf->com2 = xf->c2sav;

    /* XT-X2 midpoint CFM and the Newton system at XT, X2 */
    blmid(xf, 2);
    bldif(xf, 2);

    /* convert sensitivities wrt "T" into sensitivities wrt "1" and "2", as for the laminar part */
    for (int k = 0; k < 3; k++) {
        btrez[k] = xf->vsrez[k];
        btm[k] = xf->vsm[k]
                 + xf->vs1[0][k] * st_ms + xf->vs1[1][k] * tt_ms + xf->vs1[2][k] * dt_ms + xf->vs1[3][k] * ut_ms + xf->vs1[4][k] * xf->xt_ms;
        btr[k] = xf->vsr[k]
                 + xf->vs1[0][k] * st_re + xf->vs1[1][k] * tt_re + xf->vs1[2][k] * dt_re + xf->vs1[3][k] * ut_re + xf->vs1[4][k] * xf->xt_re;
        btx[k] = xf->vsx[k]
                 + xf->vs1[0][k] * st_xf + xf->vs1[1][k] * tt_xf + xf->vs1[2][k] * dt_xf + xf->vs1[3][k] * ut_xf + xf->vs1[4][k] * xf->xt_xf;

        bt1[k][0] = xf->vs1[0][k] * st_a1 + xf->vs1[1][k] * tt_a1 + xf->vs1[2][k] * dt_a1 + xf->vs1[3][k] * ut_a1 + xf->vs1[4][k] * xf->xt_a1;
        bt1[k][1] = xf->vs1[0][k] * st_t1 + xf->vs1[1][k] * tt_t1 + xf->vs1[2][k] * dt_t1 + xf->vs1[3][k] * ut_t1 + xf->vs1[4][k] * xf->xt_t1;
        bt1[k][2] = xf->vs1[0][k] * st_d1 + xf->vs1[1][k] * tt_d1 + xf->vs1[2][k] * dt_d1 + xf->vs1[3][k] * ut_d1 + xf->vs1[4][k] * xf->xt_d1;
        bt1[k][3] = xf->vs1[0][k] * st_u1 + xf->vs1[1][k] * tt_u1 + xf->vs1[2][k] * dt_u1 + xf->vs1[3][k] * ut_u1 + xf->vs1[4][k] * xf->xt_u1;
        bt1[k][4] = xf->vs1[0][k] * st_x1 + xf->vs1[1][k] * tt_x1 + xf->vs1[2][k] * dt_x1 + xf->vs1[3][k] * ut_x1 + xf->vs1[4][k] * xf->xt_x1;

        bt2[k][0] = xf->vs2[0][k];
        bt2[k][1] = xf->vs2[1][k]
                    + xf->vs1[0][k] * st_t2 + xf->vs1[1][k] * tt_t2 + xf->vs1[2][k] * dt_t2 + xf->vs1[3][k] * ut_t2 + xf->vs1[4][k] * xf->xt_t2;
        bt2[k][2] = xf->vs2[2][k]
                    + xf->vs1[0][k] * st_d2 + xf->vs1[1][k] * tt_d2 + xf->vs1[2][k] * dt_d2 + xf->vs1[3][k] * ut_d2 + xf->vs1[4][k] * xf->xt_d2;
        bt2[k][3] = xf->vs2[3][k]
                    + xf->vs1[0][k] * st_u2 + xf->vs1[1][k] * tt_u2 + xf->vs1[2][k] * dt_u2 + xf->vs1[3][k] * ut_u2 + xf->vs1[4][k] * xf->xt_u2;
        bt2[k][4] = xf->vs2[4][k]
                    + xf->vs1[0][k] * st_x2 + xf->vs1[1][k] * tt_x2 + xf->vs1[2][k] * dt_x2 + xf->vs1[3][k] * ut_x2 + xf->vs1[4][k] * xf->xt_x2;
    }

    /* add laminar and turbulent parts: final system in "1" and "2" variables */
    xf->vsrez[0] = btrez[0];
    xf->vsrez[1] = blrez[1] + btrez[1];
    xf->vsrez[2] = blrez[2] + btrez[2];
    xf->vsm[0] = btm[0];
    xf->vsm[1] = blm[1] + btm[1];
    xf->vsm[2] = blm[2] + btm[2];
    xf->vsr[0] = btr[0];
    xf->vsr[1] = blr[1] + btr[1];
    xf->vsr[2] = blr[2] + btr[2];
    xf->vsx[0] = btx[0];
    xf->vsx[1] = blx[1] + btx[1];
    xf->vsx[2] = blx[2] + btx[2];
    for (int l = 0; l < 5; l++) {
        xf->vs1[l][0] = bt1[0][l];
        xf->vs2[l][0] = bt2[0][l];
        xf->vs1[l][1] = bl1[1][l] + bt1[1][l];
        xf->vs2[l][1] = bl2[1][l] + bt2[1][l];
        xf->vs1[l][2] = bl1[2][l] + bt1[2][l];
        xf->vs2[l][2] = bl2[2][l] + bt2[2][l];
    }

    /* restore the "1" quantities clobbered above (the "2" ones were already restored) */
    xf->com1 = xf->c1sav;
}

void blsys(xfoil_t *xf)
{
    /* secondary BL variables and their sensitivities */
    if (xf->wake) {
        blvar(xf, 3);
        blmid(xf, 3);
    } else if (xf->turb || xf->tran) {
        blvar(xf, 2);
        blmid(xf, 2);
    } else {
        blvar(xf, 1);
        blmid(xf, 1);
    }

    /* similarity station: "1" and "2" variables are the same */
    if (xf->simi) {
        xf->com1 = xf->com2;
    }

    /* finite difference system for the current interval */
    if (xf->tran)
        trdif(xf);
    else if (xf->simi)
        bldif(xf, 0);
    else if (!xf->turb)
        bldif(xf, 1);
    else if (xf->wake)
        bldif(xf, 3);
    else if (xf->turb)
        bldif(xf, 2);

    if (xf->simi) {
        /* at the similarity station, "1" variables are really "2" variables */
        for (int k = 0; k < 4; k++) {
            for (int l = 0; l < 5; l++) {
                xf->vs2[l][k] = xf->vs1[l][k] + xf->vs2[l][k];
                xf->vs1[l][k] = 0.0;
            }
        }
    }

    /* change the system over into incompressible Uei and Mach */
    for (int k = 0; k < 4; k++) {
        /* residual derivatives wrt compressible Uec */
        const double res_u1 = xf->vs1[3][k];
        const double res_u2 = xf->vs2[3][k];
        const double res_ms = xf->vsm[k];

        /* combine with the derivatives of compressible U1, U2 = Uec(Uei, M) */
        xf->vs1[3][k] = res_u1 * xf->com1.u_uei;
        xf->vs2[3][k] = res_u2 * xf->com2.u_uei;
        xf->vsm[k] = res_u1 * xf->com1.u_ms + res_u2 * xf->com2.u_ms + res_ms;
    }
}
