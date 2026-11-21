/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  String utilities, ported from m_userio.f90.
 *
 *  Only `strip` is ported: it is the one m_userio routine reachable from the
 *  Python API that is kept (tools/fortran_callgraph.py). The console prompts
 *  aski/askr are deliberately not ported: their two reachable call sites
 *  (naca with designation <= 0, viscal with 0 iterations) get non-interactive
 *  behaviour instead (see CLAUDE.md). The other routines are unreachable.
 *
 *  Strings follow Fortran CHARACTER(len): fixed length, blank-padded, NOT
 *  NUL-terminated; the length is passed explicitly.
 */
#ifndef XFOIL_C_XUSERIO_H
#define XFOIL_C_XUSERIO_H

#include "xfoil_c.h"

/* Removes leading blanks from string[0..len-1], pads the tail with blanks, and
 * sets *ns to the number of characters up to the last non-blank one. */
XFOIL_API void strip(char *string, int len, int *ns);

#endif /* XFOIL_C_XUSERIO_H */
