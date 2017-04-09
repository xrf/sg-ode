#ifndef G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#define G_U239L9IU4I9YYJGH4M9CD6ADMGT66
#ifdef __cplusplus
extern "C" {
#endif

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
         int *restrict iwork,
         int maxnum);

#ifdef __cplusplus
}
#endif
#endif
