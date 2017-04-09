#ifndef G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#define G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#ifdef __cplusplus
extern "C" {
#endif

// TODO: remove this alias
typedef void (*U_fp)(void *, double, const double *, double *);

int ode(U_fp f, void *f_ctx, int *neqn, double *y, double *t,
        double *tout, double *relerr, double *abserr, int *
        iflag, double *work, int *iwork);

#ifdef __cplusplus
}
#endif
#endif
