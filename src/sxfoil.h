/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of s_xfoil.f90 (routines split out of xfoil to avoid circular module use):
 *  Mach/Re(CL) dependence, compressibility parameters, Cp, CL/CM integration,
 *  drag, TE panel strengths. Work on the i_xfoil globals where noted.
 *  Exported (XFOIL_API) so the tests can compare them with the Fortran originals.
 */
#ifndef XFOIL_C_SXFOIL_H
#define XFOIL_C_SXFOIL_H

#include "xfoil_c.h"
#include "xstate.h"

/* Actual minf, reinf from the unit-CL values and cls per matyp/retyp, with dM/dCL, dRe/dCL. */
XFOIL_API void mrcl(xfoil_t *xf, double cls, double *m_cls, double *r_cls);

/* Karman-Tsien parameter tklam (and tkl_msq), sonic cpstar and qstar from minf. */
XFOIL_API void comset(xfoil_t *xf);

/* Compressible cp[0..n-1] from the speeds q. */
XFOIL_API void cpcalc(int n, const double *q, double qinf, double minf, double *cp);

/* CL, CM (about xref, yref), CDp and dCL/dalpha, dCL/dM^2 from the surface vorticity. */
XFOIL_API void clcalc(int n, const double *x, const double *y, const double *gam, const double *gam_a,
                      double alfa, double minf, double qinf, double xref, double yref,
                      double *cl, double *cm, double *cdp, double *cl_alf, double *cl_msq);

/* Total drag cd (Squire-Young) and friction drag cdf. */
XFOIL_API void cdcalc(xfoil_t *xf);

/* TE gap areas ante, aste, dste, sharp, and TE panel strengths sigte, gamte (+ _a). */
XFOIL_API void tecalc(xfoil_t *xf);

#endif /* XFOIL_C_SXFOIL_H */
