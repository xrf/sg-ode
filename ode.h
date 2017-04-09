#ifndef G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#define G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif

struct Iwork {
    unsigned ns;
    bool nornd;
    unsigned k, kold;
    int isnold;
};

void ode(void (*f)(void *restrict f_ctx,
                   double t,
                   const double *restrict y,
                   double *restrict yp),
         void *restrict f_ctx,
         size_t neqn,
         double *restrict y,
         double *restrict t,
         double tout,
         double *restrict relerr,
         double *restrict abserr,
         int *restrict iflag,
         double *restrict work,
         struct Iwork *const iwork,
         int maxnum);

#ifdef __cplusplus
}
#endif
#endif
