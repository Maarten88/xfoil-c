/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Boundary-layer system, ported from m_xblsys.f90 (API-reachable routines;
 *  added incrementally). Works on the i_xbl / i_blpar globals.
 *  Closure routines: inputs by value, outputs (value and derivatives) by pointer.
 *  Exported (XFOIL_API) so the tests can compare them with the Fortran originals.
 */
#ifndef XFOIL_C_XBLSYS_H
#define XFOIL_C_XBLSYS_H

#include "xfoil_c.h"
#include "xstate.h"

/* Kinematic shape parameter hk(h, msq) (Whitfield) and derivatives. */
XFOIL_API void hkin(double h, double msq, double *hk, double *hk_h, double *hk_msq);

/* Laminar dissipation function 2 CD/H* (Falkner-Skan), and its wake version. */
XFOIL_API void dil(double hk, double rt, double *di, double *di_hk, double *di_rt);
XFOIL_API void dilw(double hk, double rt, double *di, double *di_hk, double *di_rt);

/* Laminar and turbulent H* correlations. */
XFOIL_API void hsl(double hk, double rt, double msq, double *hs, double *hs_hk, double *hs_rt, double *hs_msq);
XFOIL_API void hst(double hk, double rt, double msq, double *hs, double *hs_hk, double *hs_rt, double *hs_msq);

/* Laminar (Falkner-Skan) and turbulent (Coles; uses cffac) skin friction. */
XFOIL_API void cfl(double hk, double rt, double msq, double *cf, double *cf_hk, double *cf_rt, double *cf_msq);
XFOIL_API void cft(xfoil_t *xf, double hk, double rt, double msq, double *cf, double *cf_hk, double *cf_rt, double *cf_msq);

/* Density shape parameter (Whitfield). */
XFOIL_API void hct(double hk, double msq, double *hc, double *hc_hk, double *hc_msq);

/* Envelope e^n amplification rate dN/dx (Drela & Giles 1987, March 1991 version),
 * and the Nov 1996 version with the Orr-Sommerfeld max-amplification blend for hk > 3.5. */
XFOIL_API void dampl(double hk, double th, double rt, double *ax, double *ax_hk, double *ax_th, double *ax_rt);
XFOIL_API void dampl2(double hk, double th, double rt, double *ax, double *ax_hk, double *ax_th, double *ax_rt);

/* RMS-averaged amplification rate over interval 1..2 (dampl if idampv == 0, else
 * dampl2), plus a small term keeping dN/dx > 0 near acrit; with derivatives. */
XFOIL_API void axset(double hk1, double t1, double rt1, double a1, double hk2, double t2, double rt2, double a2,
                     double acrit, int idampv, double *ax, double *ax_hk1, double *ax_t1, double *ax_rt1,
                     double *ax_a1, double *ax_hk2, double *ax_t2, double *ax_rt2, double *ax_a2);

/* ---- station variables (work on com1/com2, see ixbl.h) ---- */

/* Sets the primary "2" variables (com2.x, ampl, s, t, d, dw, u and u's derivatives). */
XFOIL_API void blprv(xfoil_t *xf, double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei);

/* Turbulence-independent secondary "2" variables (m, r, h, v, hk, rt) from the primary ones. */
XFOIL_API void blkin(xfoil_t *xf);

/* All secondary "2" variables and their sensitivities; ityp 1 laminar, 2 turbulent, 3 wake. */
XFOIL_API void blvar(xfoil_t *xf, int ityp);

/* Midpoint skin friction cfm (and cfm_*) between stations 1 and 2. */
XFOIL_API void blmid(xfoil_t *xf, int ityp);

/* ---- difference equations (fill vs1, vs2, vsrez, vsm, vsr, vsx) ---- */

/* Newton-system coefficients and residuals of the interval 1..2: ityp 0 similarity
 * station, 1 laminar, 2 turbulent, 3 wake (transition intervals: trdif). */
XFOIL_API void bldif(xfoil_t *xf, int ityp);

/* "Dummy" BL system between the airfoil TE point and the first wake point. */
XFOIL_API void tesys(xfoil_t *xf, double cte, double tte, double dte);

/* ---- transition and system assembly ---- */

/* Checks for transition in the interval X1..X2 (2nd-order implicit amplification,
 * Dec 1994): sets ampl2, or xt and its sensitivities, and tran/trfree/trforc.
 * Returns 0 (Fortran .false.) if NaNs make convergence impossible, else 1. */
XFOIL_API int trchek2(xfoil_t *xf);
XFOIL_API int trchek(xfoil_t *xf);

/* Newton system of a transition interval: laminar X1..XT plus turbulent XT..X2. */
XFOIL_API void trdif(xfoil_t *xf);

/* Newton system of the current interval (secondary variables, then bldif/trdif),
 * converted to incompressible Uei and Mach. */
XFOIL_API void blsys(xfoil_t *xf);

#endif /* XFOIL_C_XBLSYS_H */
