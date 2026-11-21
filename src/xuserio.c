/*
 *  This file is part of xfoil_c, a C port of XFoil.
 *  Licensed under the GNU General Public License v3 or later. See LICENSE.
 *
 *  Port of strip() from m_userio.f90.
 */
#include "xuserio.h"

#include <string.h>

void strip(char *string, int len, int *ns)
{
    /* last non-blank character (-1 if none) */
    int k2 = len - 1;
    while (k2 >= 0 && string[k2] == ' ') k2--;

    /* first non-blank character; 0 if the string is blank, as the Fortran
     * zero-trip loop leaves k1 = 1 */
    int k1 = 0;
    while (k1 < k2 && string[k1] == ' ') k1++;

    *ns = k2 - k1 + 1;
    if (*ns == 0) return;

    /* shift so the first character is non-blank (Fortran overlapping
     * substring assignment = memmove), then pad the tail with blanks */
    memmove(string, string + k1, (size_t)*ns);
    memset(string + *ns, ' ', (size_t)(len - *ns));
}
