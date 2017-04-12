/* Internal header */
#ifndef G_NOEOVY8VG0HN4O8XPPLL4K2X6P1E4
#define G_NOEOVY8VG0HN4O8XPPLL4K2X6P1E4
#include <stdio.h>
#include <stdlib.h>
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

#ifdef __cplusplus
}
#endif
#endif