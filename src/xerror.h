/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Fatal errors. Fortran `stop 'msg'` becomes xfoil_stop("msg"). Inside a public API
 *  call (XFOIL_GUARDED) it records the message and jumps back to that call, which then
 *  returns normally; Python reads the message with xfoil_error() and raises XFoilError.
 *  The library state may be partly updated after an error. Outside a guard (internal
 *  routines called directly, as the tests do) it behaves like gfortran's stop: prints
 *  "STOP msg" to stderr and exits with status 0. No resources need unwinding: the
 *  library allocates no memory dynamically.
 */
#ifndef XFOIL_C_XERROR_H
#define XFOIL_C_XERROR_H

#include <setjmp.h>

void xfoil_stop(const char *msg)
#if defined(__GNUC__)
    __attribute__((noreturn))
#endif
    ;

/* Guard plumbing for XFOIL_GUARDED (not for direct use). */
int xfoil_guard_active(void);
void xfoil_guard_set(jmp_buf *jb);
void xfoil_error_clear(void);

/* Runs `stmt` (a call) so that xfoil_stop() returns here. Only the outermost public call
 * installs the guard, so API functions calling each other do not catch errors halfway. */
#define XFOIL_GUARDED(stmt)                         \
    do {                                            \
        if (xfoil_guard_active()) {                 \
            stmt;                                   \
        } else {                                    \
            jmp_buf xfoil_jb_;                      \
            xfoil_error_clear();                    \
            xfoil_guard_set(&xfoil_jb_);            \
            if (setjmp(xfoil_jb_) == 0) {           \
                stmt;                               \
            }                                       \
            xfoil_guard_set(NULL);                  \
        }                                           \
    } while (0)

#endif /* XFOIL_C_XERROR_H */
