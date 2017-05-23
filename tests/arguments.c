#include <assert.h>
#include <sg_ode.h>

#define NEQN 2

static int my_ctx;

static void f(void *ctx, double t, const double *y, double *yp)
{
    (void)t;
    assert(ctx == &my_ctx);
    assert(y);
    assert(yp);
    yp[0] =  y[1];
    yp[1] = -y[0];
}

int main(void)
{
    double abserr = 1e-5;
    double relerr = abserr;
    double bad_abserr = -abserr;
    double bad_relerr = -relerr;
    int iwork[5] = {0};
    int bad_iwork[] = {999, 999, 999, 999, 999};
    double work[100 + 21 * NEQN], y[NEQN] = {1.0, 0.0}, t = 0.0;
    const double t_out = 1.0;

    /* null arguments */
    assert(sg_ode(NULL, NULL, NEQN, y, &t, t_out,
                  relerr, abserr, 0, work, iwork) == SG_ODE_EINVAL);
    assert(sg_ode(NULL, &f, NEQN, y, NULL, t_out,
                  relerr, abserr, 0, work, iwork) == SG_ODE_EINVAL);
    assert(sg_ode(NULL, &f, NEQN, y, &t, t_out,
                  relerr, abserr, 0, NULL, iwork) == SG_ODE_EINVAL);
    assert(sg_ode(NULL, &f, NEQN, y, &t, t_out,
                  relerr, abserr, 0, work, NULL) == SG_ODE_EINVAL);

    /* invalid arguments */
    assert(sg_ode(NULL, &f, NEQN, y, &t, t_out,
                  bad_relerr, abserr, 0, work, iwork) == SG_ODE_EINVAL);
    assert(sg_ode(NULL, &f, NEQN, y, &t, t_out,
                  relerr, bad_abserr, 0, work, iwork) == SG_ODE_EINVAL);
    assert(sg_ode(NULL, &f, NEQN, y, &t, t_out,
                  relerr, abserr, 999, work, iwork) == SG_ODE_EINVAL);
    assert(sg_ode(NULL, &f, NEQN, y, &t, t_out,
                  relerr, abserr, 0, work, bad_iwork) == SG_ODE_EINVAL);

    /* nothing to do */
    assert(sg_ode(NULL, &f, NEQN, y, &t, t,
                  relerr, abserr, 0, work, iwork) == 0);

    /* normal operation */
    assert(sg_ode(&my_ctx, &f, NEQN, y, &t, t_out,
                  relerr, abserr, 0, work, iwork) == 0);
}
