/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Fatal-error handling, see xerror.h.
 */
#include "xerror.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xfoil_c.h"

/* per thread: an API call runs in the calling thread, which also reads the error afterwards */
static _Thread_local jmp_buf *active_guard = NULL;
static _Thread_local int error_pending = 0;
static _Thread_local char error_message[256];

int xfoil_guard_active(void)
{
    return active_guard != NULL;
}

void xfoil_guard_set(jmp_buf *jb)
{
    active_guard = jb;
}

void xfoil_error_clear(void)
{
    error_pending = 0;
    error_message[0] = '\0';
}

void xfoil_stop(const char *msg)
{
    fflush(stdout);
    if (active_guard) {
        snprintf(error_message, sizeof error_message, "%s", msg ? msg : "STOP");
        error_pending = 1;
        longjmp(*active_guard, 1);
    }
    if (msg) fprintf(stderr, "STOP %s\n", msg);   /* a bare Fortran `stop` prints nothing */
    exit(0);
}

int xfoil_error(char *buf, int len)
{
    if (!error_pending) return 0;
    if (buf && len > 0) snprintf(buf, (size_t)len, "%s", error_message);
    xfoil_error_clear();
    return 1;
}
