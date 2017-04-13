/* Internal header */
#ifndef G_NOEOVY8VG0HN4O8XPPLL4K2X6P1E4
#define G_NOEOVY8VG0HN4O8XPPLL4K2X6P1E4
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "inline_begin.h"
#ifdef __cplusplus
extern "C" {
#endif

static inline void *check_oom(void *ptr)
{
    if (!ptr) {
        fprintf(stderr, "Out of memory\n");
        fflush(stderr);
        abort();
    }
    return ptr;
}

static inline double min(double x, double y)
{
    return x < y ? x : y;
}

static inline double max(double x, double y)
{
    return x >= y ? x : y;
}

#if !(__STDC_VERSION__ >= 199901L || _ISOC99_SOURCE || \
      _POSIX_C_SOURCE >= 200112L ||_DEFAULT_SOURCE || \
      _BSD_SOURCE || _SVID_SOURCE)
static inline double copysign(double x, double y)
{
    return y >= 0.0 ? fabs(x) : -fabs(x);
}
#endif

#ifdef __cplusplus
}
#endif
#include "inline_end.h"
#endif
