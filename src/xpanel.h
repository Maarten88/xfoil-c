/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Panel method, ported from m_xpanel.f90 (API-reachable routines; added
 *  incrementally). Works on the i_xfoil globals (ixfoil.h).
 *  Conventions: 0-based nodes: airfoil 0..n-1, wake n..n+nw-1; a node index
 *  of -1 (Fortran 0) means "not a panel node". Logical flags are int.
 *  Exported (XFOIL_API) so the tests can compare them with the Fortran originals.
 */
#ifndef XFOIL_C_XPANEL_H
#define XFOIL_C_XPANEL_H

#include "xfoil_c.h"
#include "xstate.h"

/* Panel angles apanel[0..n-1] of the airfoil (incl. the TE panel) from x, y, nx, ny. */
XFOIL_API void apcalc(xfoil_t *xf);

/* Unit normal vectors xn, yn at the nodes of the splined curve x(s), y(s);
 * averaged at corners (equal successive s). */
XFOIL_API void ncalc(const double *x, const double *y, const double *s, int n, double *xn, double *yn);

/* Streamfunction *psi and its normal derivative *psi_ni at point (xi, yi) with
 * normal (nxi, nyi), which is node i (or -1), due to the freestream and the
 * airfoil vorticity (and sources if siglin, geometric sensitivities if geolin).
 * Sets the sensitivity globals dzdg, dzdn, dqdg, dzdm, dqdm, z_*, qtan1/2,
 * sigte, gamte, cosa, sina. */
XFOIL_API void psilin(xfoil_t *xf, int i, double xi, double yi, double nxi, double nyi,
                      double *psi, double *psi_ni, int geolin, int siglin);

/* Streamfunction and normal derivative at node i (or -1) due to the wake sources;
 * sets dzdm, dqdm for the wake nodes, and cosa, sina. */
XFOIL_API void pswlin(xfoil_t *xf, int i, double xi, double yi, double nxi, double nyi,
                      double *psi, double *psi_ni);

/* Unit vorticity distributions gamu (alpha = 0, 90 deg) from the LU-factored
 * panel system aij/aijpiv (also sets bij, qinvu, psio; lqaij = lgamu = 1). */
XFOIL_API void ggcalc(xfoil_t *xf);

/* Inviscid tangential velocity qinvu (alpha = 0, 90 deg) on the wake nodes. */
XFOIL_API void qwcalc(xfoil_t *xf);

/* Source influence matrix dij for airfoil and wake (uses bij, cij; lwdij = 1,
 * and ladij = 1 after computing the airfoil part). */
XFOIL_API void qdcalc(xfoil_t *xf);

/* Wake geometry: nw nodes n..n+nw-1 (x, y, s, nx, ny, apanel) following the
 * streamlines from the TE; sets xte, yte, lwake = 1, awake, lwdij = 0. */
XFOIL_API void xywake(xfoil_t *xf);

/* Stagnation point: panel node ist (0-based) and arc length sst, with the
 * sensitivities sst_go, sst_gp to the gam values at ist and ist+1. */
XFOIL_API void stfind(xfoil_t *xf);

/* BL station -> panel node map ipan (0-based nodes; station 0 of each side is
 * the stagnation point) and vti, iblte (0-based), nbl (counts); lipan = 1.
 * Stops (without message) if the BL arrays are too small. */
XFOIL_API void iblpan(xfoil_t *xf);

/* BL arc lengths xssi on both sides and the wake, and the TE gap profile wgap. */
XFOIL_API void xicalc(xfoil_t *xf);

/* Inviscid edge velocity uinv, uinv_a on the BL stations from qinv, qinv_a. */
XFOIL_API void uicalc(xfoil_t *xf);

/* Panel viscous tangential velocity qvis from the BL edge velocity uedg. */
XFOIL_API void qvfue(xfoil_t *xf);

/* Inviscid panel tangential velocity qinv, qinv_a for the current alfa. */
XFOIL_API void qiset(xfoil_t *xf);

/* gam = qvis and gam_a = qinv_a on the airfoil. */
XFOIL_API void gamqv(xfoil_t *xf);

/* Moves the stagnation point to the new panel (from gam) and shifts the BL
 * variables (ctau, thet, dstr, uedg, itran) accordingly; resets mass. */
XFOIL_API void stmove(xfoil_t *xf);

/* uedg = uinv + all source influence (dij * mass). */
XFOIL_API void ueset(xfoil_t *xf);

#endif /* XFOIL_C_XPANEL_H */
