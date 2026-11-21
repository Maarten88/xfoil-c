/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Storage layout of the BL Newton system's mass-coupling blocks VM (optimization O2).
 *
 *  Exact build (default): the generated layout xf->vm[iv][jv][k] = Fortran VM(k,jv,iv), byte-identical
 *  to gfortran's (the tests compare it with the Fortran library).
 *  XFOIL_FAST build: the same memory is used as [iv][k][jv], so that a row k of a block line iv is
 *  contiguous in jv; blsolv's eliminations then run over contiguous memory and vectorise. The
 *  arithmetic is the same, so the results are the same; only the memory order of VM differs.
 */
#ifndef XFOIL_C_XLAYOUT_H
#define XFOIL_C_XLAYOUT_H

#include "xstate.h"

#if XFOIL_FAST
/* VM(k, jv, iv) of the Fortran: row k (0..2) of block line iv, column jv */
#define XF_VM(xf, k, jv, iv) (((double (*)[3][IZX])(void *)(xf)->vm)[iv][k][jv])
#else
#define XF_VM(xf, k, jv, iv) ((xf)->vm[iv][jv][k])
#endif

#endif /* XFOIL_C_XLAYOUT_H */
