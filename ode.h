#ifndef G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#define G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#include <stdbool.h>
#include <math.h>
#include "vector.h"
#ifdef __cplusplus
extern "C" {
#endif

struct Ode {
    struct SgVectorDriver drv;
    double alpha[12], beta[12], sig[13], v[12], w[12], g[13], psi[12];
    double x, h, hold, told, delsgn;
    double *yy, *wt, *p, *yp, *ypout, *phi[16];
    unsigned ns;
    // invariant: k >= 1 && k < 13
    unsigned k;
    unsigned kold;
    // "iflag_sign_old": whether the user-provided iflag was positive (controls
    // whether the solver is allowed to overshoot and interpolate)
    bool isnold;
    // probably means "no_round" and has something to do with rounding
    bool nornd;
    bool phase1;
    bool start;
};

void ode_init(struct Ode *self, struct SgVectorDriver drv);

void ode_del(struct Ode *self);

void ode(struct Ode *self,
         void (*f)(void *f_ctx,
                   double t,
                   const SgVector *restrict y,
                   SgVector *restrict yp),
         void *restrict f_ctx,
         struct SgVectorDriver drv,
         SgVector *restrict y,
         double *restrict t,
         double tout,
         double *restrict relerr,
         double *restrict abserr,
         unsigned maxnum,
         int *restrict iflag);

#ifdef __cplusplus
}
#endif
#endif
