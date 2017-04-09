#ifndef G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#define G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#ifdef __cplusplus
extern "C" {
#endif

int ode(void (*f)(void *f_ctx, double t, const double *y, double *yp),
        void *f_ctx, int neqn, double *y, double *t,
        double tout, double *relerr, double *abserr, int *
        iflag, double *work, int *iwork, int maxnum);

#ifdef __cplusplus
}
#endif
#endif
