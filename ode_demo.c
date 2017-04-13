/** @file
    Example program for the Shampine-Gordon ODE solver.

    Solves the defining equations of Jacobian elliptic functions

        y₀′ =     y₁ y₂    y₀(0) = 0
        y₁′ =    -y₀ y₂    y₁(0) = 1
        y₂′ = -k² y₀ y₁    y₂(0) = 1

    for t[i] ∈ {0, 5, 10, 15, …, 60}.  The analytic solution is

        y₀ = sn(t / k²)
        y₁ = cn(t / k²)
        y₂ = dn(t / k²)

    in this case we use k² = 0.51.
*/
#include <stdio.h>
#include "ode.h"

void f(void *ctx, double t, const double *y_vec, double *yp_vec)
{
    static const double ksq = 0.51;
    const double *const y = (const double *)y_vec;
    double *const yp = (double *)yp_vec;
    (void)ctx;
    (void)t;
    yp[0] = y[1] * y[2];
    yp[1] = -y[0] * y[2];
    yp[2] = -ksq * y[0] * y[1];
}

int main(void)
{
    static const size_t neqn = 3;
    static const double relerr = 1.0e-9;
    static const double abserr = 1.0e-16;
    int iflag = 1;
    double t;
    int i;
    int iwork[5] = {0};
    double y[3], work[100 + 21 * 3];

    printf("neqn=%zu relerr=%g abserr=%g iflag=%i\n",
           neqn, relerr, abserr, iflag);

#define dump() \
    printf("t=%2.0f y=[%17.10e %17.10e %17.10e] iflag=%i\n", \
           t, y[0], y[1], y[2], iflag);

    t = 0.0;
    y[0] = 0.0;
    y[1] = 1.0;
    y[2] = 1.0;
    dump();

    for (i = 1; i <= 12; ++i) {
        double tout = 5.0 * i;
        int r;
        r = sg_ode(NULL, &f, neqn, y, &t, tout, relerr, abserr, 0, work, iwork);
        iflag = r == 0 ? 2 : iflag;
        dump();
    }
    return 0;
}
