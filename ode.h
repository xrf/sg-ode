#ifndef G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#define G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#include <stdbool.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

struct Ode {
    // invariant: k >= 1 && k < 13
    unsigned ns, k, kold;
    int isnold; // what is this??
    double alpha[12], beta[12], sig[13], v[12], w[12], g[13], psi[12];
    double x, h, hold, told, delsgn;
    bool nornd, phase1, start;
};

#define ODE_INITIALIZER {65535, 65535, 65535, 32767, {NAN}, {NAN}, {NAN}, {NAN}, {NAN}, {NAN}, {NAN}, NAN, NAN, NAN, NAN, NAN, false, false, false}

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
         struct Ode *const iwork,
         int maxnum);

#ifdef __cplusplus
}
#endif
#endif
