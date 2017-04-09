#ifndef G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#define G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#include <stdbool.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

struct Ode {
    double alpha[12], beta[12], sig[13], v[12], w[12], g[13], psi[12];
    double x, h, hold, told, delsgn;
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
         struct Ode *const self,
         int maxnum);

#ifdef __cplusplus
}
#endif
#endif
