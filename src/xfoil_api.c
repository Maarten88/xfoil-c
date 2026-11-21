/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Public API, ported from api.f90 (bind(c) routines, all arguments by pointer).
 *  The full contract is in API_SNAPSHOT.md.
 */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xstate.h"
#include "sxfoil.h"
#include "xerror.h"
#include "xfoil.h"
#include "xfoil_c.h"
#include "xfortran.h"
#include "xgdes.h"
#include "xgeom.h"
#include "xoper.h"
#include "xspline.h"

#ifndef XFOIL_C_VERSION
#  error "XFOIL_C_VERSION must be defined by the build (CMakeLists.txt)"
#endif

const char *xfoil_version(void)
{
    return XFOIL_C_VERSION;
}

const char *xfoil_build_info(void)
{
    return XFOIL_BUILD_INFO;
}

xfoil_t *xfoil_new(void)
{
    xfoil_t *xf = calloc(1, sizeof *xf);     /* zero pages come lazily from the OS: cheap */
    if (xf) xfoil_state_init(xf);
    return xf;
}

void xfoil_free(xfoil_t *xf)
{
    free(xf);
}

void set_print(xfoil_t *xf, const bool *setting)
{
    xf->show_output = *setting ? 1 : 0;   /* Fortran logical(4) */
    xfoil_print = xf->show_output;
}

bool get_print(xfoil_t *xf)
{
    return xf->show_output != 0;
}

static void set_airfoil_impl(xfoil_t *xf, const double *x_in, const double *y_in, const int *n_in)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */
    xf->scr.flap_base_valid = 0;    /* new airfoil: no flap deflection to undo */
    xf->scr.flap_angle = 0.0;
    const double angtol = 40.0;

    /* the Fortran wrote past the buffer for more than IBX points (clean-up C6a) */
    if (*n_in > IBX) xfoil_stop("SET_AIRFOIL: too many coordinate points (more than IBX = 2560)");
    if (*n_in < 3) xfoil_stop("SET_AIRFOIL: at least 3 coordinate points are needed");

    for (int i = 0; i < *n_in; i++) {
        xf->xb[i] = x_in[i];
        xf->yb[i] = y_in[i];
    }
    xf->nb = *n_in;

    /* calculate airfoil area assuming counterclockwise ordering */
    double area = 0.0;
    for (int i = 0; i < xf->nb; i++) {
        const int ip = (i == xf->nb - 1) ? 0 : i + 1;
        area = area + 0.5 * (xf->yb[i] + xf->yb[ip]) * (xf->xb[i] - xf->xb[ip]);
    }

    if (area >= 0.0) {
        if (xfoil_print) printf("\n Number of input coordinate points:%4d\n Counterclockwise ordering\n", xf->nb);
    } else {
        /* if area is negative (clockwise order), reverse coordinate order */
        if (xfoil_print) printf("\n Number of input coordinate points:%4d\n Clockwise ordering\n", xf->nb);
        for (int i = 0; i < xf->nb / 2; i++) {
            const double xtmp = xf->xb[xf->nb - 1 - i];
            const double ytmp = xf->yb[xf->nb - 1 - i];
            xf->xb[xf->nb - 1 - i] = xf->xb[i];
            xf->yb[xf->nb - 1 - i] = xf->yb[i];
            xf->xb[i] = xtmp;
            xf->yb[i] = ytmp;
        }
    }

    if (xf->lnorm) {
        norm(xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb);
        if (xfoil_print) printf("\n Airfoil has been normalized\n");
    }

    scalc(xf->xb, xf->yb, xf->sb, xf->nb);
    segspl(xf->xb, xf->xbp, xf->sb, xf->nb);
    segspl(xf->yb, xf->ybp, xf->sb, xf->nb);

    geopar(xf->xb, xf->xbp, xf->yb, xf->ybp, xf->sb, xf->nb, xf->w1, &xf->sble, &xf->chordb, &xf->areab, &xf->radble, &xf->angbte, &xf->ei11ba, &xf->ei22ba,
           &xf->apx1ba, &xf->apx2ba, &xf->ei11bt, &xf->ei22bt, &xf->apx1bt, &xf->apx2bt, &xf->thickb, &xf->cambrb);

    double xble = seval(xf->sble, xf->xb, xf->xbp, xf->sb, xf->nb);
    double yble = seval(xf->sble, xf->yb, xf->ybp, xf->sb, xf->nb);
    const double xbte = 0.5 * (xf->xb[0] + xf->xb[xf->nb - 1]);
    const double ybte = 0.5 * (xf->yb[0] + xf->yb[xf->nb - 1]);

    if (xfoil_print) {
        printf("\n  LE  x,y  =%10.5f%10.5f  |   Chord =%10.5f\n", xble, yble, xf->chordb);
        printf("  TE  x,y  =%10.5f%10.5f  |\n", xbte, ybte);
    }

    /* wipe out old flap hinge location */
    xf->xbf = 0.0;
    xf->ybf = 0.0;
    xf->lbflap = 0;

    abcopy(xf, 1);

    int imax;
    double amax;
    cang(xf->x, xf->y, xf->n, 0, &imax, &amax);
    if (fabs(amax) > angtol && xfoil_print) {
        printf("\n WARNING: Poor input coordinate distribution\n");
        printf("          Excessive panel angle%7.1f  at i =%4d\n", amax, imax + 1);
        printf("          Repaneling with PANE and/or PPAR suggested\n");
        printf("           (doing GDES,CADD before repaneling _may_\n");
        printf("            improve excessively coarse LE spacing\n");
    }
}

int get_n_coords(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    return xf->nb;
}

void get_airfoil(xfoil_t *xf, double *x_out, double *y_out, const int *n_out)
{
    xfoil_print = xf->show_output;
    for (int i = 0; i < *n_out; i++) {
        x_out[i] = xf->xb[i];
        y_out[i] = xf->yb[i];
    }
}

static void set_naca_impl(xfoil_t *xf, const int *spec)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */
    xf->scr.flap_base_valid = 0;    /* new airfoil: no flap deflection to undo */
    xf->scr.flap_angle = 0.0;
    if (*spec <= 0) {
        /* Fortran writes to unit 0 (stderr) regardless of show_output */
        fprintf(stderr, " Invalid NACA specifier. Specify a NACA 4 or 5 series airfoil code.\n");
    } else {
        naca(xf, *spec);
    }
}

void set_reynolds(xfoil_t *xf, const double *re)
{
    xfoil_print = xf->show_output;
    xf->lvisc = (*re == 0.0) ? 0 : 1;
    xf->reinf1 = *re;
    xf->lblini = 0;
    xf->lipan = 0;
    xf->lvconv = 0;
}

double get_reynolds(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    return xf->reinf1;
}

static void set_mach_impl(xfoil_t *xf, const double *m)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */
    if (*m != xf->minf) {
        xf->minf = *m;
        xf->minf1 = *m;

        comset(xf);
        if (xf->minf > 0.0 && xfoil_print)
            printf("\n Sonic Cp =%10.2f      Sonic Q/Qinf =%10.3f\n\n", xf->cpstar, xf->qstar / xf->qinf);

        cpcalc(xf->n, xf->qinv, xf->qinf, xf->minf, xf->cpi);
        if (xf->lvisc) cpcalc(xf->n + xf->nw, xf->qvis, xf->qinf, xf->minf, xf->cpv);
        clcalc(xf->n, xf->x, xf->y, xf->gam, xf->gam_a, xf->alfa, xf->minf, xf->qinf, xf->xcmref, xf->ycmref, &xf->cl, &xf->cm, &xf->cdp, &xf->cl_alf, &xf->cl_msq);
        cdcalc(xf);
    }
    xf->lvconv = 0;
}

double get_mach(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    return xf->minf;
}

void set_xtr(xfoil_t *xf, const double *xtr_top, const double *xtr_bot)
{
    xfoil_print = xf->show_output;
    xf->xstrip[0] = *xtr_top;
    xf->xstrip[1] = *xtr_bot;
    xf->lvconv = 0;
}

void get_xtr(xfoil_t *xf, double *xtr_top, double *xtr_bot)
{
    xfoil_print = xf->show_output;
    *xtr_top = xf->xstrip[0];
    *xtr_bot = xf->xstrip[1];
}

void set_n_crit(xfoil_t *xf, const double *n_crit)
{
    xfoil_print = xf->show_output;
    xf->acrit[0] = *n_crit;
    xf->acrit[1] = *n_crit;
    xf->lvconv = 0;
}

double get_n_crit(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    return xf->acrit[0];
}

void set_max_iter(xfoil_t *xf, const int *max_iter)
{
    xfoil_print = xf->show_output;
    xf->itmax = *max_iter;
}

int get_max_iter(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    return xf->itmax;
}

void reset_bls(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    xf->lblini = 0;
    xf->lipan = 0;
    if (xfoil_print) printf(" BLs will be initialized on next point\n");
}

static void repanel_impl(xfoil_t *xf, const int *n_panel, const double *cv_par, const double *cte_ratio, const double *ctr_ratio,
             const double *xs_ref1, const double *xs_ref2, const double *xp_ref1, const double *xp_ref2)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */
    xf->npan = (*n_panel < IQX - 6) ? *n_panel : IQX - 6;   /* integer min */
    xf->cvpar = *cv_par;
    xf->cterat = *cte_ratio;
    xf->ctrrat = *ctr_ratio;
    xf->xsref1 = *xs_ref1;
    xf->xsref2 = *xs_ref2;
    xf->xpref1 = *xp_ref1;
    xf->xpref2 = *xp_ref2;

    pangen(xf, 1);
    if (xf->n > 0) {
        int imax;
        double amax;
        cang(xf->x, xf->y, xf->n, 1, &imax, &amax);
    }
}


/* ---- analysis ---------------------------------------------------------- */

#define P2(v) __builtin_powi((v), 2)

static void alfa_impl(xfoil_t *xf, const double *a_input, double *cl_out, double *cd_out, double *cm_out, double *cp_out,
          double *xtrtop_out, double *xtrbot_out, bool *conv)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */
    xf->adeg = *a_input;

    xf->alfa = *a_input * xf->dtor;
    xf->qinf = 1.0;
    xf->lalfa = 1;
    specal(xf);
    if (fabs(xf->alfa - xf->awake) > 1.0e-5) xf->lwake = 0;
    if (fabs(xf->alfa - xf->avisc) > 1.0e-5) xf->lvconv = 0;
    if (fabs(xf->minf - xf->mvisc) > 1.0e-5) xf->lvconv = 0;

    if (xf->lvisc) {
        *conv = viscal(xf, xf->itmax) != 0;
        *conv = xf->lvconv && *conv;
    } else {
        *conv = true;
    }

    *cl_out = xf->cl;
    *cd_out = xf->cd;
    *cm_out = xf->cm;

    fcpmin(xf);
    *cp_out = xf->cpmn;

    *xtrtop_out = xf->xoctr[0];
    *xtrbot_out = xf->xoctr[1];
    xf->scr.cp_valid = *conv;
}

static void cl_impl(xfoil_t *xf, const double *cl_input, double *a_out, double *cd_out, double *cm_out, double *cp_out,
        double *xtrtop_out, double *xtrbot_out, bool *conv)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */
    xf->clspec = *cl_input;
    xf->alfa = 0.0;
    xf->qinf = 1.0;
    xf->lalfa = 0;
    speccl(xf);
    xf->adeg = xf->alfa / xf->dtor;
    if (fabs(xf->alfa - xf->awake) > 1.0e-5) xf->lwake = 0;
    if (fabs(xf->alfa - xf->avisc) > 1.0e-5) xf->lvconv = 0;
    if (fabs(xf->minf - xf->mvisc) > 1.0e-5) xf->lvconv = 0;

    if (xf->lvisc) {
        *conv = viscal(xf, xf->itmax) != 0;
        *conv = xf->lvconv && *conv;
    } else {
        *conv = true;
    }

    *a_out = xf->alfa / xf->dtor;
    *cd_out = xf->cd;
    *cm_out = xf->cm;

    fcpmin(xf);
    *cp_out = xf->cpmn;

    *xtrtop_out = xf->xoctr[0];
    *xtrbot_out = xf->xoctr[1];
    xf->scr.cp_valid = *conv;
}

/* Store one sequence point and its convergence flag (shared tail of aseq/cseq). conv_arr[i]
 * keeps the caller's value when the point is inviscid, as in Fortran (it is then overwritten
 * with .true. anyway). */
static void store_seq_point(xfoil_t *xf, int i, double *a_arr, double *cl_arr, double *cd_arr, double *cm_arr,
                            double *cp_arr, bool *conv_arr)
{
    xf->adeg = xf->alfa / xf->dtor;

    a_arr[i] = xf->adeg;
    cl_arr[i] = xf->cl;
    cd_arr[i] = xf->cd;
    cm_arr[i] = xf->cm;

    fcpmin(xf);
    cp_arr[i] = xf->cpmn;

    if ((xf->lvconv && conv_arr[i]) || !xf->lvisc) {
        conv_arr[i] = true;
    } else if (xf->lvisc && !(xf->lvconv && conv_arr[i])) {
        conv_arr[i] = false;
    }
}

static void aseq_impl(xfoil_t *xf, const double *a_start, const double *a_end, const int *n_step, double *a_arr, double *cl_arr,
          double *cd_arr, double *cm_arr, double *cp_arr, bool *conv_arr)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */
    const double a0 = *a_start * xf->dtor;
    /* divides by n_step (not n_step - 1): a_end itself is not reached, as in Fortran */
    const double da = (*a_end - *a_start) / (double)*n_step * xf->dtor;

    xf->lalfa = 1;

    for (int i = 0; i < *n_step; i++) {
        xf->alfa = a0 + da * (double)i;
        if (fabs(xf->alfa - xf->awake) > 1.0e-5) xf->lwake = 0;
        if (fabs(xf->alfa - xf->avisc) > 1.0e-5) xf->lvconv = 0;
        if (fabs(xf->minf - xf->mvisc) > 1.0e-5) xf->lvconv = 0;
        specal(xf);
        const int itmaxs = xf->itmax + 5;
        if (xf->lvisc) conv_arr[i] = viscal(xf, itmaxs) != 0;
        store_seq_point(xf, i, a_arr, cl_arr, cd_arr, cm_arr, cp_arr, conv_arr);
    }
    if (*n_step > 0) xf->scr.cp_valid = conv_arr[*n_step - 1];
}

static void cseq_impl(xfoil_t *xf, const double *cl_start, const double *cl_end, const int *n_step, double *a_arr, double *cl_arr,
          double *cd_arr, double *cm_arr, double *cp_arr, bool *conv_arr)
{
    xf->scr.cp_valid = 0;     /* get_cp: NaN until a converged analysis (clean-up C6c) */
    const double cl0 = *cl_start;
    const double dcl = (*cl_end - *cl_start) / (double)*n_step;
    xf->lalfa = 0;

    for (int i = 0; i < *n_step; i++) {
        xf->clspec = cl0 + dcl * (double)i;
        speccl(xf);
        if (fabs(xf->alfa - xf->awake) > 1.0e-5) xf->lwake = 0;
        if (fabs(xf->alfa - xf->avisc) > 1.0e-5) xf->lvconv = 0;
        if (fabs(xf->minf - xf->mvisc) > 1.0e-5) xf->lvconv = 0;

        const int itmaxs = xf->itmax + 5;
        if (xf->lvisc) conv_arr[i] = viscal(xf, itmaxs) != 0;
        store_seq_point(xf, i, a_arr, cl_arr, cd_arr, cm_arr, cp_arr, conv_arr);
    }
    if (*n_step > 0) xf->scr.cp_valid = conv_arr[*n_step - 1];
}

int get_n_cp(xfoil_t *xf)
{
    xfoil_print = xf->show_output;
    return xf->n;
}

static void get_cp_impl(xfoil_t *xf, double *x_out, double *y_out, double *cp_out, const int *n_points)
{
    (void)n_points;   /* Fortran ignores it and writes N values (Phase 12) */
    comset(xf);
    const double beta = sqrt(1.0 - P2(xf->minf));
    const double bfac = 0.5 * P2(xf->minf) / (1.0 + beta);

    /* NaN unless the last analysed point converged and neither the geometry nor the Mach number
     * changed since (Fortran returned the stale or unconverged values; clean-up C6c) */
    for (int i = 0; i < xf->n; i++) {
        const double cpinc = 1.0 - P2(xf->gam[i] / xf->qinf);
        const double den = beta + bfac * cpinc;
        cp_out[i] = xf->scr.cp_valid ? cpinc / den : NAN;
        x_out[i] = xf->x[i];
        y_out[i] = xf->y[i];
    }
}

/* ---- public entry points: run inside the error guard (xerror.h) ---- */

void set_airfoil(xfoil_t *xf, const double *x_in, const double *y_in, const int *n_in)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(set_airfoil_impl(xf, x_in, y_in, n_in));
}

void set_naca(xfoil_t *xf, const int *spec)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(set_naca_impl(xf, spec));
}

void set_mach(xfoil_t *xf, const double *m)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(set_mach_impl(xf, m));
}

void repanel(xfoil_t *xf, const int *n_panel, const double *cv_par, const double *cte_ratio, const double *ctr_ratio,
             const double *xs_ref1, const double *xs_ref2, const double *xp_ref1, const double *xp_ref2)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(repanel_impl(xf, n_panel, cv_par, cte_ratio, ctr_ratio, xs_ref1, xs_ref2, xp_ref1, xp_ref2));
}

void alfa(xfoil_t *xf, const double *a_input, double *cl_out, double *cd_out, double *cm_out, double *cp_out,
          double *xtrtop_out, double *xtrbot_out, bool *conv)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(alfa_impl(xf, a_input, cl_out, cd_out, cm_out, cp_out, xtrtop_out, xtrbot_out, conv));
}

void cl(xfoil_t *xf, const double *cl_input, double *a_out, double *cd_out, double *cm_out, double *cp_out,
        double *xtrtop_out, double *xtrbot_out, bool *conv)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(cl_impl(xf, cl_input, a_out, cd_out, cm_out, cp_out, xtrtop_out, xtrbot_out, conv));
}

void aseq(xfoil_t *xf, const double *a_start, const double *a_end, const int *n_step, double *a_arr, double *cl_arr,
          double *cd_arr, double *cm_arr, double *cp_arr, bool *conv_arr)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(aseq_impl(xf, a_start, a_end, n_step, a_arr, cl_arr, cd_arr, cm_arr, cp_arr, conv_arr));
}

void cseq(xfoil_t *xf, const double *cl_start, const double *cl_end, const int *n_step, double *a_arr, double *cl_arr,
          double *cd_arr, double *cm_arr, double *cp_arr, bool *conv_arr)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(cseq_impl(xf, cl_start, cl_end, n_step, a_arr, cl_arr, cd_arr, cm_arr, cp_arr, conv_arr));
}

void get_cp(xfoil_t *xf, double *x_out, double *y_out, double *cp_out, const int *n_points)
{
    xfoil_print = xf->show_output;
    XFOIL_GUARDED(get_cp_impl(xf, x_out, y_out, cp_out, n_points));
}
