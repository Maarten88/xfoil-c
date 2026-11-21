/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Per-instance state that is not a Fortran module variable: work arrays of single routines
 *  that were static (Fortran local arrays live in static memory there, so they keep their
 *  contents between calls; clean-up C5d) and API bookkeeping (C6c). Part of the context, so
 *  every XFoil instance and thread has its own; zero at creation like statics.
 */
#ifndef XFOIL_C_XSCRATCH_H
#define XFOIL_C_XSCRATCH_H

#include "ixfoil.h"

typedef struct xfoil_scratch {
    /* setbl (xbl.c): sensitivities of the "1"/"2" station variables wrt the system lines, saved Ue */
    double d1_m[2 * IVX], d2_m[2 * IVX], u1_m[2 * IVX], u2_m[2 * IVX];
    double ule1_m[2 * IVX], ule2_m[2 * IVX], ute1_m[2 * IVX], ute2_m[2 * IVX];
    double usav[2][IVX];
    /* update (xbl.c): new Ue / panel speeds and their sensitivities wrt alpha or CL */
    double unew[2][IVX], u_ac[2][IVX], qnew[IQX], q_ac[IQX];
    /* get_geometry (xfoil_ext.c): geopar work array */
    double geo_t[IZX];
    /* flap (xfoil_ext.c): the undeflected buffer airfoil and the current deflection, so that
     * flap(angle) sets the angle instead of adding to it */
    double flap_base_xb[IBX], flap_base_yb[IBX];
    int flap_base_nb, flap_base_valid;
    double flap_angle;
    /* API: the last analysed point converged and nothing changed the geometry or Mach since;
     * otherwise get_cp returns NaN for cp (clean-up C6c) */
    int cp_valid;
} xfoil_scratch_t;

#endif /* XFOIL_C_XSCRATCH_H */
